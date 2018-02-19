/**
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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "LagrangeLambdaConstraint.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<LagrangeLambdaConstraint>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for Lagrange multiplier problem.");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetic vector");
  params.addCoupledVar("mag_z", 0.0, "The z component of the magnetic vector");
  params.addRequiredCoupledVar("lambda", "The lagrange multiplier");
  params.addParam<Real>("eps",1.0,"eps");
  return params;
}

LagrangeLambdaConstraint::LagrangeLambdaConstraint(const InputParameters & parameters)
  :Kernel(parameters),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _lambda(coupledValue("lambda")),
  _eps(getParam<Real>("eps"))
{
}

Real
LagrangeLambdaConstraint::computeQpResidual()
{
  return _test[_i][_qp] * ((_mag_x[_qp] * _mag_x[_qp] + _mag_y[_qp] * _mag_y[_qp] + _mag_z[_qp] * _mag_z[_qp] - 1.0 - _eps * _lambda[_qp]));
}

Real
LagrangeLambdaConstraint::computeQpJacobian()
{
  return - _eps * _test[_i][_qp];
}

Real
LagrangeLambdaConstraint::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _mag_x_var)
  {
    return 2.0 * _test[_i][_qp] * _phi[_j][_qp] * _mag_x[_qp];
  }
  else if (jvar == _mag_y_var)
  {
    return 2.0 * _test[_i][_qp] * _phi[_j][_qp] * _mag_y[_qp];
  }
  else if (jvar == _mag_z_var)
  {
    return 2.0 * _test[_i][_qp] * _phi[_j][_qp] * _mag_y[_qp];
  }
  else
    return 0.0;
}
