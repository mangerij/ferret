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

#include "BandGapAuxTiO2.h"

registerMooseObject("FerretApp", BandGapAuxTiO2);

template <>
InputParameters validParams<BandGapAuxTiO2>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculates the changes to local band gap due to the elastic stress fields.");
  params.addParam<Real>("relaxed_energy", 0.0,"relaxed energy");
  params.addParam<Real>("biaxial_stress_rate", 0.0, "biaxial stress rate");
  params.addParam<Real>("uniaxial_stress_rate", 0.0, "uniaxial stress rate");
  return params;
}


BandGapAuxTiO2::BandGapAuxTiO2(const InputParameters & parameters) :
  AuxKernel(parameters),
   _stress(getMaterialProperty<RankTwoTensor>("stress")),
   _ba(getParam<Real>("biaxial_stress_rate")),
   _bc(getParam<Real>("uniaxial_stress_rate")),
   _E0(getParam<Real>("relaxed_energy"))
{
}

Real
BandGapAuxTiO2::computeValue()

{
    return _E0 + _ba*_stress[_qp](0,0)+_bc*_stress[_qp](2,2);
}


