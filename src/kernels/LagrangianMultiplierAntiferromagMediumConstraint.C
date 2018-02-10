/**
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANT_lambda[_qp]ILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "LagrangianMultiplierAntiferromagMediumConstraint.h"

template<>
InputParameters validParams<LagrangianMultiplierAntiferromagMediumConstraint>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution to enforce M$*$M =1. Note: untested");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the vari_lambda[_qp]le this kernel acts in. (0 for x, 1 for y, 2 for z, 3 for lambda)");
  params.addRequiredCoupledVar("antiferromag_L_x", "The x component of the antiferromagnetic vector");
  params.addRequiredCoupledVar("antiferromag_L_y", "The y component of the antiferromagnetic vector");
  params.addRequiredCoupledVar("antiferromag_L_z", "The z component of the antiferromagnetic vector");
  params.addRequiredCoupledVar("lambda", "The Lagrange multiplier field vari_lambda[_qp]le");
  params.addRequiredParam<Real>("epsilon", "The epsilon parameter to enforce nonzero diagonals.");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

LagrangianMultiplierAntiferromagMediumConstraint::LagrangianMultiplierAntiferromagMediumConstraint(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _antiferromag_L_x_var(coupled("antiferromag_L_x")),
  _antiferromag_L_y_var(coupled("antiferromag_L_y")),
  _antiferromag_L_z_var(coupled("antiferromag_L_z")),
  _antiferromag_L_x(coupledValue("antiferromag_L_x")),
  _antiferromag_L_y(coupledValue("antiferromag_L_y")),
  _antiferromag_L_z(coupledValue("antiferromag_L_z")),
  _lambda_var(coupled("lambda")),
  _lambda(coupledValue("lambda")),
  _epsilon(getParam<Real>("epsilon")),
  _len_scale(getParam<Real>("len_scale"))
{
}


Real
LagrangianMultiplierAntiferromagMediumConstraint::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (-2*_lambda[_qp]*_antiferromag_L_x[_qp] - 2*std::pow(_lambda[_qp],2)*_antiferromag_L_x[_qp] + 2*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_x[_qp],3) + 2*std::pow(_lambda[_qp],2)*_antiferromag_L_x[_qp]*std::pow(_antiferromag_L_y[_qp],2) + 2*std::pow(_lambda[_qp],2)*_antiferromag_L_x[_qp]*std::pow(_antiferromag_L_z[_qp],2));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (-2*_lambda[_qp]*_antiferromag_L_y[_qp] - 2*std::pow(_lambda[_qp],2)*_antiferromag_L_y[_qp] + 2*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_x[_qp],2)*_antiferromag_L_y[_qp] + 2*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_y[_qp],3) + 2*std::pow(_lambda[_qp],2)*_antiferromag_L_y[_qp]*std::pow(_antiferromag_L_z[_qp],2));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (-2*_lambda[_qp]*_antiferromag_L_z[_qp] - 2*std::pow(_lambda[_qp],2)*_antiferromag_L_z[_qp] + 2*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_x[_qp],2)*_antiferromag_L_z[_qp] + 2*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_y[_qp],2)*_antiferromag_L_z[_qp] + 2*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_z[_qp],3));
  }
  else if (_component == 3)
  {
    return _test[_i][_qp] * (1 + _lambda[_qp] - std::pow(_antiferromag_L_x[_qp],2) - 2*_lambda[_qp]*std::pow(_antiferromag_L_x[_qp],2) + _lambda[_qp]*std::pow(_antiferromag_L_x[_qp],4) - std::pow(_antiferromag_L_y[_qp],2) - 2*_lambda[_qp]*std::pow(_antiferromag_L_y[_qp],2) + 2*_lambda[_qp]*std::pow(_antiferromag_L_x[_qp],2)*std::pow(_antiferromag_L_y[_qp],2) + _lambda[_qp]*std::pow(_antiferromag_L_y[_qp],4) - std::pow(_antiferromag_L_z[_qp],2) - 2*_lambda[_qp]*std::pow(_antiferromag_L_z[_qp],2) + 2*_lambda[_qp]*std::pow(_antiferromag_L_x[_qp],2)*std::pow(_antiferromag_L_z[_qp],2) + 2*_lambda[_qp]*std::pow(_antiferromag_L_y[_qp],2)*std::pow(_antiferromag_L_z[_qp],2) + _lambda[_qp]*std::pow(_antiferromag_L_z[_qp],4));
  }
  else
  {
    return 0.0;
  }
}

Real
LagrangianMultiplierAntiferromagMediumConstraint::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (-2*_lambda[_qp] - 2*std::pow(_lambda[_qp],2) + 6*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_x[_qp],2) + 2*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_y[_qp],2) + 2*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_z[_qp],2));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (-2*_lambda[_qp] - 2*std::pow(_lambda[_qp],2) + 2*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_x[_qp],2) + 6*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_y[_qp],2) + 2*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_z[_qp],2));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (-2*_lambda[_qp] - 2*std::pow(_lambda[_qp],2) + 2*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_x[_qp],2) + 2*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_y[_qp],2) + 6*std::pow(_lambda[_qp],2)*std::pow(_antiferromag_L_z[_qp],2));
  }
  else if (_component == 3)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (1 - 2*std::pow(_antiferromag_L_x[_qp],2) + std::pow(_antiferromag_L_x[_qp],4) - 2*std::pow(_antiferromag_L_y[_qp],2) + 2*std::pow(_antiferromag_L_x[_qp],2)*std::pow(_antiferromag_L_y[_qp],2) + std::pow(_antiferromag_L_y[_qp],4) - 2*std::pow(_antiferromag_L_z[_qp],2) + 2*std::pow(_antiferromag_L_x[_qp],2)*std::pow(_antiferromag_L_z[_qp],2) + 2*std::pow(_antiferromag_L_y[_qp],2)*std::pow(_antiferromag_L_z[_qp],2) + std::pow(_antiferromag_L_z[_qp],4));
  }
  else
  {
    return 0.0;
  }
}

Real
LagrangianMultiplierAntiferromagMediumConstraint::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiferromag_L_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*std::pow(_lambda[_qp],2)*_antiferromag_L_x[_qp]*_antiferromag_L_y[_qp]);
    }
    else if (jvar == _antiferromag_L_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*std::pow(_lambda[_qp],2)*_antiferromag_L_x[_qp]*_antiferromag_L_z[_qp]);
    }
    else if (jvar == _lambda_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-2*_antiferromag_L_x[_qp] - 4*_lambda[_qp]*_antiferromag_L_x[_qp] + 4*_lambda[_qp]*std::pow(_antiferromag_L_x[_qp],3) + 4*_lambda[_qp]*_antiferromag_L_x[_qp]*std::pow(_antiferromag_L_y[_qp],2) + 4*_lambda[_qp]*_antiferromag_L_x[_qp]*std::pow(_antiferromag_L_z[_qp],2));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _antiferromag_L_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*std::pow(_lambda[_qp],2)*_antiferromag_L_x[_qp]*_antiferromag_L_y[_qp]);
    }
    else if (jvar == _antiferromag_L_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*std::pow(_lambda[_qp],2)*_antiferromag_L_y[_qp]*_antiferromag_L_z[_qp]);
    }
    else if (jvar == _lambda_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-2*_antiferromag_L_y[_qp] - 4*_lambda[_qp]*_antiferromag_L_y[_qp] + 4*_lambda[_qp]*std::pow(_antiferromag_L_x[_qp],2)*_antiferromag_L_y[_qp] + 4*_lambda[_qp]*std::pow(_antiferromag_L_y[_qp],3) + 4*_lambda[_qp]*_antiferromag_L_y[_qp]*std::pow(_antiferromag_L_z[_qp],2));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _antiferromag_L_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*std::pow(_lambda[_qp],2)*_antiferromag_L_x[_qp]*_antiferromag_L_z[_qp]);
    }
    else if (jvar == _antiferromag_L_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (4*std::pow(_lambda[_qp],2)*_antiferromag_L_y[_qp]*_antiferromag_L_z[_qp]);
    }
    else if (jvar == _lambda_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-2*_antiferromag_L_z[_qp] - 4*_lambda[_qp]*_antiferromag_L_z[_qp] + 4*_lambda[_qp]*std::pow(_antiferromag_L_x[_qp],2)*_antiferromag_L_z[_qp] + 4*_lambda[_qp]*std::pow(_antiferromag_L_y[_qp],2)*_antiferromag_L_z[_qp] + 4*_lambda[_qp]*std::pow(_antiferromag_L_z[_qp],3));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 3)
  {
    if (jvar == _antiferromag_L_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-2*_antiferromag_L_x[_qp] - 4*_lambda[_qp]*_antiferromag_L_x[_qp] + 4*_lambda[_qp]*std::pow(_antiferromag_L_x[_qp],3) + 4*_lambda[_qp]*_antiferromag_L_x[_qp]*std::pow(_antiferromag_L_y[_qp],2) + 4*_lambda[_qp]*_antiferromag_L_x[_qp]*std::pow(_antiferromag_L_z[_qp],2));
    }
    else if (jvar == _antiferromag_L_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-2*_antiferromag_L_y[_qp] - 4*_lambda[_qp]*_antiferromag_L_y[_qp] + 4*_lambda[_qp]*std::pow(_antiferromag_L_x[_qp],2)*_antiferromag_L_y[_qp] + 4*_lambda[_qp]*std::pow(_antiferromag_L_y[_qp],3) + 4*_lambda[_qp]*_antiferromag_L_y[_qp]*std::pow(_antiferromag_L_z[_qp],2));
    }
    else if (jvar == _antiferromag_L_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (-2*_antiferromag_L_z[_qp] - 4*_lambda[_qp]*_antiferromag_L_z[_qp] + 4*_lambda[_qp]*std::pow(_antiferromag_L_x[_qp],2)*_antiferromag_L_z[_qp] + 4*_lambda[_qp]*std::pow(_antiferromag_L_y[_qp],2)*_antiferromag_L_z[_qp] + 4*_lambda[_qp]*std::pow(_antiferromag_L_z[_qp],3));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
