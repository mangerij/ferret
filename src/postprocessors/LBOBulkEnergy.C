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

#include "LBOBulkEnergy.h"

template<>
InputParameters validParams<LBOBulkEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over local bulk energy density for LBO-type pseudo-uniaxial ferroelectrics.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha2", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha3", "The coefficients of the Landau expansion");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

LBOBulkEnergy::LBOBulkEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
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
LBOBulkEnergy::computeQpIntegral()
{
  return (- 0.5 * _alpha1 * _polar_z[_qp] * _polar_z[_qp] + 0.25 * _alpha2 * _polar_z[_qp] * _polar_z[_qp] * _polar_z[_qp] * _polar_z[_qp] + 0.5 * _alpha3 * (_polar_x[_qp] * _polar_x[_qp] + _polar_y[_qp] * _polar_y[_qp])) * std::pow(_len_scale,3);
}
