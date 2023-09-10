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

#include "AFMTotalEnergyDensity.h"
registerMooseObject("FerretApp", AFMTotalEnergyDensity);

InputParameters AFMTotalEnergyDensity::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the sum of energy densities");
  params.addRequiredCoupledVar("Edmi", "The DMI energy density");
  params.addRequiredCoupledVar("Esupexch", "The superexchange energy density");
  params.addRequiredCoupledVar("Enlexch", "The non-local exchange stiffness energy density");
  params.addRequiredCoupledVar("Eepa1", "The easy-plane anisotropy energy density for sublattice 1");
  params.addRequiredCoupledVar("Eepa2", "The easy-plane anisotropy energy density for sublattice 2");
  params.addRequiredCoupledVar("Eca1", "The cubic anisotropy energy density for sublattice 1");
  params.addRequiredCoupledVar("Eca2", "The cubic anisotropy energy density for sublattice 2");
  return params;
}


AFMTotalEnergyDensity::AFMTotalEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
   _Edmi(coupledValue("Edmi")),
   _Esupexch(coupledValue("Esupexch")),
   _Enlexch(coupledValue("Enlexch")),
   _Eepa1(coupledValue("Eepa1")),
   _Eepa2(coupledValue("Eepa2")),
   _Eca1(coupledValue("Eca1")),
   _Eca2(coupledValue("Eca2"))
{
}

Real
AFMTotalEnergyDensity::computeValue()
{
  return _Edmi[_qp] + _Esupexch[_qp] + _Enlexch[_qp] + _Eepa1[_qp] + _Eepa2[_qp] + _Eca1[_qp] + _Eca2[_qp];
}
