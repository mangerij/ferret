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

#include "ChangeInRefractiveIndexElectro.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

template<>

InputParameters validParams<ChangeInRefractiveIndexElectro>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("index_i", "first index of the beta vector");
  params.addRequiredParam<unsigned int>("index_j", "second index of the beta vector");
  params.addRequiredParam<unsigned int>("index_k", "first index of the delta_beta vector");
  params.addRequiredParam<unsigned int>("index_l", "second index of the delta_beta vector");
  return params;
}


ChangeInRefractiveIndexElectro::ChangeInRefractiveIndexElectro(const InputParameters & parameters) :
  AuxKernel(parameters),
   _index_i(getParam<unsigned int>("index_i")),
   _index_j(getParam<unsigned int>("index_j")),
   _index_k(getParam<unsigned int>("index_k")),
   _index_l(getParam<unsigned int>("index_l")),
   _indicatrix(getMaterialProperty<RankTwoTensor>("indicatrix")),
   _delta_indicatrix_electro(getMaterialProperty<RankTwoTensor>("delta_indicatrix_electro")) //or need RealTensorValue
{
}

Real
ChangeInRefractiveIndexElectro::computeValue()
{
  // the diagonals are related to the B1, B2, B3 terms in rotated indicatrix
//std::pow(  (1.0 / ( _beta_tensor[_qp](_index_i, _index_j)  ) ), 3.0)
  return - 0.5 * std::pow(1.0 / std::pow(_indicatrix[_qp](_index_i, _index_j), 0.5) , 3.0) *  _delta_indicatrix_electro[_qp](_index_k, _index_l);
}
