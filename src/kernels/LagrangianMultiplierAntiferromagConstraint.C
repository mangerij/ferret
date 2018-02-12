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

#include "LagrangianMultiplierAntiferromagConstraint.h"

template<>
InputParameters validParams<LagrangianMultiplierAntiferromagConstraint>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution to enforce Lagrange mulitplier-based constraints of M$*$M = 1");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z, 3 for lambda)");
  params.addRequiredCoupledVar("mag_x", "The x component of the antiferromagnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the antiferromagnetic vector");
  params.addCoupledVar("mag_z", 0.0, "The z component of the antiferromagnetic vector");
  params.addRequiredCoupledVar("lambda", "The Lagrange multiplier field variable");
  params.addRequiredParam<Real>("epsilon", "The epsilon parameter to enforce nonzero diagonals.");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

LagrangianMultiplierAntiferromagConstraint::LagrangianMultiplierAntiferromagConstraint(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _lambda_var(coupled("lambda")),
  _lambda(coupledValue("lambda")),
  _epsilon(getParam<Real>("epsilon")),
  _len_scale(getParam<Real>("len_scale"))
{
}


Real
LagrangianMultiplierAntiferromagConstraint::computeQpResidual()
{
  if (_component == 0)
  {
    return -_test[_i][_qp] * (2.0*_lambda[_qp]*_mag_x[_qp]);
  }
  else if (_component == 1)
  {
    return -_test[_i][_qp] * (2.0*_lambda[_qp]*_mag_y[_qp]);
  }
  else if (_component == 2)
  {
    return -_test[_i][_qp] * (2.0*_lambda[_qp]*_mag_z[_qp]);
  }
  else if (_component == 3)
  {
    return _test[_i][_qp] * (1.0 - (std::pow(_mag_x[_qp],2) + std::pow(_mag_y[_qp],2) + std::pow(_mag_z[_qp],2)));
  }
  else
  {
    return 0.0;
  }
}

Real
LagrangianMultiplierAntiferromagConstraint::computeQpJacobian()
{
  if (_component == 0)
  {
    return -_test[_i][_qp] * _phi[_j][_qp] * (2.0*_lambda[_qp]);
  }
  else if (_component == 1)
  {
    return -_test[_i][_qp] * _phi[_j][_qp] * (2.0*_lambda[_qp]);
  }
  else if (_component == 2)
  {
    return -_test[_i][_qp] * _phi[_j][_qp] * (2.0*_lambda[_qp]);
  }
  else if (_component == 3)
  {
    return 0.0;
  }
  else
  {
    return 0.0;
  }
}

Real
LagrangianMultiplierAntiferromagConstraint::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _lambda_var)
    {
      return -_test[_i][_qp] * _phi[_j][_qp] * (2.0 * _mag_x[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _lambda_var)
    {
      return -_test[_i][_qp] * _phi[_j][_qp] * (2.0 * _mag_y[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _lambda_var)
    {
      return -_test[_i][_qp] * _phi[_j][_qp] * (2.0 * _mag_z[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 3)
  {
    if (jvar == _mag_x_var)
    {
      return -2.0 * _test[_i][_qp] * _phi[_j][_qp] * _mag_x[_qp];
    }
    else if (jvar == _mag_y_var)
    {
      return -2.0 * _test[_i][_qp] * _phi[_j][_qp] * _mag_y[_qp];
    }
    else if (jvar == _mag_z_var)
    {
      return -2.0 * _test[_i][_qp] * _phi[_j][_qp] * _mag_z[_qp];
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
