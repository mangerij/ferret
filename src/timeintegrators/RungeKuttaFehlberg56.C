/*
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "RungeKuttaFehlberg56.h"
#include "NonlinearSystemBase.h"
#include "FEProblem.h"
#include "PetscSupport.h"

registerMooseObject("MooseApp", RungeKuttaFehlberg56);

defineLegacyParams(RungeKuttaFehlberg56);

InputParameters
RungeKuttaFehlberg56::validParams()
{
  InputParameters params = TimeIntegrator::validParams();
  params.addClassDescription(
      "Sixth-order Runge Kutta method with embedded fifth order solution");
  return params;
}

// Initialize static data
const Real RungeKuttaFehlberg56::_c[RungeKuttaFehlberg56::_n_stages] = {0.0,1.0/6.0,4.0/15.0,2.0/3.0,4.0/5.0,1.0,0.0,1.0};

const Real RungeKuttaFehlberg56::_a[RungeKuttaFehlberg56::_n_stages][RungeKuttaFehlberg56::_n_stages] = 
    {
      {1.0,      0.0,         0.0, 0.0, 0.0, 0.0, 0.0, 0.0},                                 //A1i
      {1.0/6.0,  0.0,         0.0, 0.0, 0.0, 0.0, 0.0, 0.0},                                 //A2i
      {4.0/75.0, 16.0/75.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0},                                 //A3i
      {5.0/6.0,  -16.0/6.0, 15.0/6.0, 0.0, 0.0, 0.0, 0.0, 0.0},                              //A4i
      {-40.0/25.0, 144.0/25.0, -100.0/25.0, 16.0/25.0, 0.0, 0.0, 0.0, 0.0},                  //A5i
      {722.0/640.0, -2304.0/640.0, 2035.0/640.0, -88.0/640.0, 275.0/640.0, 0.0, 0.0, 0.0},   //A6i
      {-22.0/1280.0, 0.0, 55.0/1280.0, -88.0/1280.0, 55.0/1280.0, 0.0, 0.0, 0.0},           //A7i
      {186.0/1280.0, -4608.0/1280.0, 4015.0/1280.0, -88.0/1280.0, 495.0/1280.0, 0.0, 0.0, 0.0}         //A8i
    };

RungeKuttaFehlberg56::RungeKuttaFehlberg56(const InputParameters & parameters)
  : TimeIntegrator(parameters), _stage(1)
{
  mooseInfo("RungeKuttaFehlberg56 and other multistage TimeIntegrators are known not to work with "
            "Materials/AuxKernels that accumulate 'state' and should be used with caution.");

  // Name the stage residuals "residual_stage1", "residual_stage2", etc.
  for (unsigned int stage = 0; stage < _n_stages; ++stage)
  {
    std::ostringstream oss;
    oss << "residual_stage" << stage + 1;
    _stage_residuals[stage] = &(_nl.addVector(oss.str(), false, GHOSTED));
  }
}

void
RungeKuttaFehlberg56::computeTimeDerivatives()
{
  // We are multiplying by the method coefficients in postResidual(), so
  // the time derivatives are of the same form at every stage although
  // the current solution varies depending on the stage.
  if (!_sys.solutionUDot())
    mooseError("RungeKuttaFehlberg56: Time derivative of solution (`u_dot`) is not stored. Please set "
               "uDotRequested() to true in FEProblemBase befor requesting `u_dot`.");

  NumericVector<Number> & u_dot = *_sys.solutionUDot();
  u_dot = *_solution;
  computeTimeDerivativeHelper(u_dot, _solution_old);
  u_dot.close();
  _du_dot_du = 1. / _dt;
}

void
RungeKuttaFehlberg56::computeADTimeDerivatives(DualReal & ad_u_dot,
                                       const dof_id_type & dof,
                                       DualReal & /*ad_u_dotdot*/) const
{
  computeTimeDerivativeHelper(ad_u_dot, _solution_old(dof));
}

void
RungeKuttaFehlberg56::solve()
{
  // Time at end of step
  Real time_old = _fe_problem.timeOld();

  // Reset iteration counts
  _n_nonlinear_iterations = 0;
  _n_linear_iterations = 0;

  // A for-loop would increment _stage too far, so we use an extra
  // loop counter.
  for (unsigned int current_stage = 1; current_stage <= _n_stages; ++current_stage)
  {
    // Set the current stage value
    _stage = current_stage;

    // This ensures that all the Output objects in the OutputWarehouse
    // have had solveSetup() called, and sets the default solver
    // parameters for PETSc.
    _fe_problem.initPetscOutput();

    _console << "Stage " << _stage << std::endl;

    // Set the time for this stage
    _fe_problem.time() = time_old + _c[_stage - 1] * _dt;

    // Do the solve
    _fe_problem.getNonlinearSystemBase().system().solve();

    // Update the iteration counts
    _n_nonlinear_iterations += getNumNonlinearIterationsLastSolve();
    _n_linear_iterations += getNumLinearIterationsLastSolve();

    // Abort time step immediately on stage failure - see TimeIntegrator doc page
    if (!_fe_problem.converged())
      return;
  }
}

void
RungeKuttaFehlberg56::postResidual(NumericVector<Number> & residual)
{
  // Error if _stage got messed up somehow.
  if (_stage > _n_stages)
    // the explicit cast prevents strange compiler weirdness with the static
    // const variable and the variadic mooseError function
    mooseError("RungeKuttaFehlberg56::postResidual(): Member variable _stage can only have values 1-",
               (unsigned int)_n_stages,
               ".");

  // In the standard RK notation, the residual of stage 1 of s is given by:
  //
  // R := M*(Y_i - y_n)/dt - \sum_{j=1}^s a_{ij} * f(t_n + c_j*dt, Y_j) = 0
  //
  // where:
  // .) M is the mass matrix
  // .) Y_i is the stage solution
  // .) dt is the timestep, and is accounted for in the _Re_time residual.
  // .) f are the "non-time" residuals evaluated for a given stage solution.
  // .) The minus signs are already "baked in" to the residuals and so do not appear below.

  // Store this stage's non-time residual.  We are calling operator=
  // here, and that calls close().
  *_stage_residuals[_stage - 1] = _Re_non_time;

  // Build up the residual for this stage.
  residual.add(1., _Re_time);
  for (unsigned int j = 0; j < _stage; ++j)
    residual.add(_a[_stage - 1][j], *_stage_residuals[j]);
  residual.close();
}
