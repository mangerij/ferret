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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ThermoelectricZTAux.h"

registerMooseObject("FerretApp", ThermoelectricZTAux);

InputParameters ThermoelectricZTAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates thermoelectric figure of merit");
  params.addCoupledVar("T", "Temperature variable in Kelvin");

  return params;
}

ThermoelectricZTAux::ThermoelectricZTAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _T_var(coupled("T")),
  _T(coupledValue("T")),
  _T_grad(coupledGradient("T")),
  _ecC(getMaterialProperty<Real>("ecC")),
  _sbC(getMaterialProperty<Real>("sbC")),
  _thC(getMaterialProperty<Real>("thC"))
{
}

Real
ThermoelectricZTAux::computeValue()
{
    return _ecC[_qp] * _sbC[_qp] * _sbC[_qp]  / _thC[_qp] * _T[_qp];
}
