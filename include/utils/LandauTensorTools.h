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

#ifndef LANDAUTENSORTOOLS_H
#define LANDAUTENSORTOOLS_H

template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
template <typename>
class RankFourTensorTempl;
typedef RankFourTensorTempl<Real> RankFourTensor;
//class RankSixTensor;
//class RankEightTensor;

namespace LandauTensorTools
{
  Real landauTwoProduct(const RankTwoTensor & aij, const RealVectorValue & p);
  Real landauFourProduct(const RankFourTensor & aijkl, const RealVectorValue & p);
 // Real landauSixProduct(const RankSixTensor & aijklmn, const RealVectorValue & p);
  Real landauTwoProductDerivative(const RankTwoTensor & aij, unsigned int k, const RealVectorValue & p);
  Real landauFourProductDerivative(const RankFourTensor & aijkl, unsigned int m, const RealVectorValue & p);
//  Real landauSixProductDerivative(const RankSixTensor & aijklmn, unsigned int m, const RealVectorValue & p);
}

#endif //LANDAUTENSORTOOLS_H
