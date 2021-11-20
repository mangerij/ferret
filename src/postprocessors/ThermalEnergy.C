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

#include "ThermalEnergy.h"

registerMooseObject("FerretApp", ThermalEnergy);

InputParameters ThermalEnergy::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral whose integrand is 3$/$2$*$k$*$T where k is a constant (Boltzmann perhaps).");
  params.addRequiredCoupledVar("temperature", "The local temperature");
  params.addParam<Real>("const", 1.0, "thermal relation constant (e.g. Boltzmann constant)");
  return params;
}

ThermalEnergy::ThermalEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _temperature(coupledValue("temperature")),
   _const(getParam<Real>("const"))
{
}

Real
ThermalEnergy::computeQpIntegral()
{
  return (3.0/2.0) * _const * _temperature[_qp];
}
