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

#include "LBOBulkEnergyDeriv.h"

class LBOBulkEnergyDeriv;

registerMooseObject("FerretApp", LBOBulkEnergyDeriv);

template<>
InputParameters validParams<LBOBulkEnergyDeriv>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to local bulk energy of LBO, where B can be Ti, Ta, etc.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha2", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha3", "The coefficients of the Landau expansion");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

LBOBulkEnergyDeriv::LBOBulkEnergyDeriv(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _alpha1(getParam<Real>("alpha1")),
   _alpha2(getParam<Real>("alpha2")),
   _alpha3(getParam<Real>("alpha3")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
LBOBulkEnergyDeriv::computeQpResidual()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  if (_component == 0)
  {
    return _test[_i][_qp] * ( _alpha3 * w(0) );
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * ( _alpha3 * w(1) );
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (- _alpha1  * w(2) +  _alpha2  * w(2)*w(2)*w(2) );
  }
  else
    return 0.0;
}

Real
LBOBulkEnergyDeriv::computeQpJacobian()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  if (_component == 0)
  {
    return _test[_i][_qp] * _alpha3 * _phi[_j][_qp];
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _alpha3 * _phi[_j][_qp];
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (- _alpha1 + 3 * _alpha2  * w(2) * w(2) ) * _phi[_j][_qp];
  }
  else
    return 0.0;
}
