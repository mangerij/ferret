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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "TimeLambdaConstraintLLG.h"
#include "MooseVariable.h"

registerMooseObject("FerretApp", TimeLambdaConstraintLLG);

template<>
InputParameters validParams<TimeLambdaConstraintLLG>()
{
  InputParameters params = validParams<TimeKernel>();
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  params.addRequiredCoupledVar("lambda", "The Lagrange multiplier");
  params.addRequiredParam<Real>("eps", "eps");
  return params;
}

TimeLambdaConstraintLLG::TimeLambdaConstraintLLG(const InputParameters & parameters) :
    TimeKernel(parameters),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _lambda_var(coupled("lambda")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _lambda(coupledValue("lambda")),
  _mag_x_dot(coupledDot("mag_x")),
  _mag_y_dot(coupledDot("mag_y")),
  _mag_z_dot(coupledDot("mag_z")),
  _mag_x_d_dot(coupledDotDu("mag_x")),
  _mag_y_d_dot(coupledDotDu("mag_y")),
  _mag_z_d_dot(coupledDotDu("mag_z")),
  _eps(getParam<Real>("eps"))
{
}

Real
TimeLambdaConstraintLLG::computeQpResidual()
{
  return _test[_i][_qp]*(-_eps*_lambda[_qp]+_mag_x[_qp]*_mag_x_dot[_qp]+_mag_y[_qp]*_mag_y_dot[_qp]+_mag_z[_qp]*_mag_z_dot[_qp]);
}

Real
TimeLambdaConstraintLLG::computeQpJacobian()
{
  return -_test[_i][_qp]*_eps*_phi[_j][_qp];
}

Real
TimeLambdaConstraintLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _mag_x_var)
  {
    return _test[_i][_qp]*_phi[_j][_qp]*(_mag_x_dot[_qp]+_mag_x[_qp]*_mag_x_d_dot[_qp]);
  }
  else if (jvar == _mag_y_var)
  {
    return _test[_i][_qp]*_phi[_j][_qp]*(_mag_y_dot[_qp]+_mag_y[_qp]*_mag_y_d_dot[_qp]);
  }
  else if (jvar == _mag_z_var)
  {
    return _test[_i][_qp]*_phi[_j][_qp]*(_mag_z_dot[_qp]+_mag_z[_qp]*_mag_z_d_dot[_qp]);
  }
  else
  {
    return 0.0;
  }
}

