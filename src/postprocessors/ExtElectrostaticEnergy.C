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

#include "ExtElectrostaticEnergy.h"

template<>
InputParameters validParams<ExtElectrostaticEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("potential_E_int", "The internal electric potential");
  params.addRequiredParam<Real>("permittivity", "permittivity");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

ExtElectrostaticEnergy::ExtElectrostaticEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _potential_E_int_grad(coupledGradient("potential_E_int")),
   _permittivity(getParam<Real>("permittivity")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
ExtElectrostaticEnergy::computeQpIntegral()
{
  return 0.5 * _permittivity * _potential_E_int_grad[_qp] * _potential_E_int_grad[_qp] * std::pow(_len_scale, 2.0);
}
