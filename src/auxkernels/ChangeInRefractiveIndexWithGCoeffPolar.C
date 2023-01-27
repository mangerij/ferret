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

#include "ChangeInRefractiveIndexWithGCoeffPolar.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ChangeInRefractiveIndexWithGCoeffPolar);

InputParameters ChangeInRefractiveIndexWithGCoeffPolar::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the changes to local refractive index due to the polar-optic effect.");
  params.addRequiredParam<unsigned int>("index_i", "first index of the beta vector");
  params.addRequiredParam<unsigned int>("index_j", "second index of the beta vector");
  params.addRequiredParam<unsigned int>("index_k", "first index of the delta_beta vector");
  params.addRequiredParam<unsigned int>("index_l", "second index of the delta_beta vector");
  return params;
}


ChangeInRefractiveIndexWithGCoeffPolar::ChangeInRefractiveIndexWithGCoeffPolar(const InputParameters & parameters) :
  AuxKernel(parameters),
   _index_i(getParam<unsigned int>("index_i")),
   _index_j(getParam<unsigned int>("index_j")),
   _index_k(getParam<unsigned int>("index_k")),
   _index_l(getParam<unsigned int>("index_l")),
   _indicatrix(getMaterialProperty<RankTwoTensor>("indicatrix")),
   _delta_indicatrix(getMaterialProperty<RankTwoTensor>("delta_indicatrix")),
   _delta_gPO_tensor(getMaterialProperty<RankTwoTensor>("delta_gPO_tensor"))
{
}

Real
ChangeInRefractiveIndexWithGCoeffPolar::computeValue()
{
  // the diagonals are related to the B1, B2, B3 terms in rotated indicatrix
  // Likely we will need to load in the dielectric constant representation of the inverse refractive index squared instead of this constant number here.
  return - 0.5 * Utility::pow<3>(1.0 / std::pow(_indicatrix[_qp](_index_i, _index_j), 0.5) ) *  (- _delta_indicatrix[_qp](_index_k, _index_l) - _delta_gPO_tensor[_qp](_index_k, _index_l));
}
