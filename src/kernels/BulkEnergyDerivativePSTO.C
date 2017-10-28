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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "BulkEnergyDerivativePSTO.h"
#include<cmath>

template<>
InputParameters validParams<BulkEnergyDerivativePSTO>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredParam<Real>("alpha1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha2", "");
  params.addRequiredParam<Real>("alpha3", "");
  params.addRequiredParam<Real>("alpha4", "");
  params.addRequiredParam<Real>("alpha5", "");
  params.addRequiredParam<Real>("x1", "The coupling terms");
  params.addRequiredParam<Real>("x2", "");
  params.addRequiredParam<Real>("x3", "");
  params.addRequiredParam<Real>("x4", "");
  params.addRequiredParam<Real>("x5", "");
  params.addRequiredParam<Real>("x6", "");
  params.addRequiredParam<Real>("epsilon", "Constant value for strain");
  params.addRequiredParam<Real>("T", "Temperature");
  return params;
}

BulkEnergyDerivativePSTO::BulkEnergyDerivativePSTO(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _alpha1(getParam<Real>("alpha1")),
   _alpha2(getParam<Real>("alpha2")),
   _alpha3(getParam<Real>("alpha3")),
   _alpha4(getParam<Real>("alpha4")),
   _alpha5(getParam<Real>("alpha5")),
   _x1(getParam<Real>("x1")),
   _x2(getParam<Real>("x2")),
   _x3(getParam<Real>("x3")),
   _x4(getParam<Real>("x4")),
   _x5(getParam<Real>("x5")),
   _x6(getParam<Real>("x6")),
   _epsilon(getParam<Real>("epsilon")),
   _T(getParam<Real>("T"))
{
}

Real
BulkEnergyDerivativePSTO::computeQpResidual()
{
  if (_component == 0)
  {
    return (_alpha1 * (_T - 120.0) * (2.0 * _polar_x[_qp]) + _alpha2 * 4.0 * std::pow(_polar_x[_qp], 3.0) + 
_alpha3 * 2.0 * _polar_x[_qp] * std::pow(_polar_y[_qp], 2.0) + 6.0 * _alpha4 * std::pow(_polar_x[_qp], 5.0) + 
_alpha5 * (4.0 * std::pow(_polar_x[_qp], 3.0) * std::pow(_polar_y[_qp], 2.0) + 2.0 * _polar_x[_qp] * std::pow(_polar_y[_qp], 4.0)) + 
(_x1 * ( 2.0 * _polar_x[_qp]) + _x2 * (4.0 * std::pow(_polar_x[_qp], 3.0)) + 
_x3 * 2.0 * _polar_x[_qp] * std::pow(_polar_y[_qp], 2.0)) * _epsilon + (_x4 * (2.0 * _polar_x[_qp]) + 
_x5 * ( 4.0 * std::pow(_polar_x[_qp], 3.0)) + _x6 * 2.0 * _polar_x[_qp] * std::pow(_polar_y[_qp], 2.0)) * std::pow(_epsilon, 2.0))  * _test[_i][_qp];
  }
  else if (_component == 1)
  {
   return (_alpha1 * (_T - 120.0) * (2.0 * _polar_y[_qp] ) + _alpha2 * 4.0 * std::pow(_polar_y[_qp], 3.0) + 
_alpha3 * 2.0 * _polar_y[_qp] * std::pow(_polar_x[_qp], 2.0) + 6.0 * _alpha4 *std::pow(_polar_y[_qp], 5.0) + 
_alpha5 * (4.0 * std::pow(_polar_y[_qp], 3.0) * std::pow(_polar_x[_qp], 2.0) + 2.0 * _polar_y[_qp] * std::pow(_polar_x[_qp], 4.0)) + 
(_x1 * ( 2.0 * _polar_y[_qp]) + _x2 * (4.0 * std::pow(_polar_y[_qp], 3.0)) + 
_x3 * 2.0 * _polar_y[_qp] * std::pow(_polar_x[_qp], 2.0)) * _epsilon + (_x4 * (2.0 * _polar_y[_qp]) + 
_x5 * ( 4.0 * std::pow(_polar_y[_qp], 3.0)) + _x6 * 2.0 * _polar_y[_qp] * std::pow(_polar_x[_qp], 2.0)) * std::pow(_epsilon, 2.0))  * _test[_i][_qp];
  }
  else 
  {
    return 0.0;
  }
}

Real
BulkEnergyDerivativePSTO::computeQpJacobian()
{
  if (_component == 0)
  {
    return (_alpha1 * (_T - 120.0) * 2.0 + _alpha2 * 12.0 * std::pow(_polar_x[_qp], 2.0) + 
_alpha3 * 2.0 * std::pow(_polar_y[_qp], 2.0) + 30.0 * _alpha4 * std::pow(_polar_x[_qp], 4.0) + 
_alpha5 * (12.0 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0) + 2.0 * std::pow(_polar_y[_qp], 4.0)) + 
(_x1 * 2.0 + _x2 * 12.0 * std::pow(_polar_x[_qp], 2.0) + 
_x3 * 2.0 * std::pow(_polar_y[_qp], 2.0)) * _epsilon + (_x4 * 2.0 + 
_x5 * 12.0 * std::pow(_polar_x[_qp], 2.0) + _x6 * 2.0 * std::pow(_polar_y[_qp], 2.0)) * std::pow(_epsilon, 2.0)) * _phi[_j][_qp] * _test[_i][_qp];    
  }
  else if (_component == 1)
  {
    return (_alpha1 * (_T - 120.0) * 2.0 + _alpha2 * 12.0 * std::pow(_polar_y[_qp], 2.0) + 
_alpha3 * 2.0 * std::pow(_polar_x[_qp], 2.0) + 30.0 * _alpha4 * std::pow(_polar_y[_qp], 4.0) + 
_alpha5 * (12.0 * std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_x[_qp], 2.0) + 2.0 * std::pow(_polar_x[_qp], 4.0)) + 
(_x1 * 2.0 + _x2 * 12.0 * std::pow(_polar_y[_qp], 2.0) + 
_x3 * 2.0 * std::pow(_polar_x[_qp], 2.0)) * _epsilon + (_x4 * 2.0 + 
_x5 * 12.0 * std::pow(_polar_y[_qp], 2.0) + _x6 * 2.0 * std::pow(_polar_x[_qp], 2.0)) * std::pow(_epsilon, 2.0)) * _phi[_j][_qp] * _test[_i][_qp];
  }
  else 
  {
    return 0.0;
  }
}

Real
BulkEnergyDerivativePSTO::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 1)
  {
    if (jvar == _polar_x_var)
    {
      return (_alpha3 * 4.0 * _polar_x[_qp] * _polar_y[_qp] + 
_alpha5 * (8.0 * std::pow(_polar_x[_qp], 3.0) * _polar_y[_qp] + 8.0 * _polar_x[_qp] * std::pow(_polar_y[_qp], 3.0)) + 
(_x3 * 4.0 * _polar_x[_qp] * _polar_y[_qp]) * _epsilon +
(_x6 * 4.0 * _polar_x[_qp] * _polar_y[_qp]) * std::pow(_epsilon, 2.0)) * _test[_i][_qp] * _phi[_j][_qp];
    }
    else
      return 0.0;
  }
  
  else if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return (_alpha3 * 4.0 * _polar_x[_qp] * _polar_y[_qp] + 
_alpha5 * (8.0 * std::pow(_polar_x[_qp], 3.0) * _polar_y[_qp] + 8.0 * _polar_x[_qp] * std::pow(_polar_y[_qp], 3.0)) + 
(_x3 * 4.0 * _polar_x[_qp] * _polar_y[_qp]) * _epsilon +
(_x6 * 4.0 * _polar_x[_qp] * _polar_y[_qp]) * std::pow(_epsilon, 2.0)) * _test[_i][_qp] * _phi[_j][_qp];
    }
    else
      return 0.0;
  }
  else 
    return 0.0;
}




