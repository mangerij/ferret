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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "HoleDensityAux.h"
registerMooseObject("FerretApp", HoleDensityAux);

InputParameters HoleDensityAux::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<Real>("Ev", "Property name of the conduction band energy (J)");
  params.addRequiredParam<Real>("Nv", "Effective DOS of the conduction band (T=298)");
  params.addRequiredParam<Real>("T", "temperature (K)");
  params.addRequiredParam<Real>("Kb", "Boltzmann Constant (aJ/K)");
  params.addRequiredParam<Real>("q", "eV (aJ)");
  params.addRequiredCoupledVar("potential_E_int","E");
  return params;
}


HoleDensityAux::HoleDensityAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _Ev(getParam<Real>("Ev")),
  _Nv(getParam<Real>("Nv")),
  _T(getParam<Real>("T")),
  _Kb(getParam<Real>("Kb")),
  _q(getParam<Real>("q")),
  _potential_E_int(coupledValue("potential_E_int"))
{
}

Real
HoleDensityAux::computeValue()
{
  return (_Nv * std::exp((( _Ev - (_q * _potential_E_int[_qp])) / (_Kb * _T))));
}
