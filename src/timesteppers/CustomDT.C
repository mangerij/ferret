/**
 * @file   CustomDT.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Fri Nov 15 13:17:57 2013
 *
 * @brief  Demonstrate the "ungluing of linear and nonlinear residual
 *
 *
 */
#include "CustomDT.h"
#include "FEProblem.h"
//libMesh includes
#include "libmesh/implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/transient_system.h"
// C++ Includes
#include <iomanip>
#include<limits>

template<>
InputParameters validParams<CustomDT>()
{
  InputParameters params = validParams<TimeStepper>();
  params.addParam<Real>("dt", 1., "The initial time step size.");
  params.addParam<Real>("increase_rate",2.0, "the rate the step size increases");
  return params;
}

CustomDT::CustomDT(const std::string & name, InputParameters parameters) :
  TimeStepper(name, parameters),
    _increase_rate(getParam<Real>("increase_rate"))
{}

Real
CustomDT::computeInitialDT()
{
  return getParam<Real>("dt");
}

Real
CustomDT::computeDT()
{
  Real new_dt;
  if(_t_step==3) new_dt=getCurrentDT()*_increase_rate;
  else new_dt=getCurrentDT();
  return new_dt;
}
