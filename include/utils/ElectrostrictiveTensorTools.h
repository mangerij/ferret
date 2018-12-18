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

#ifndef ELECTROSTRICTIVETENSORTOOLS_H
#define ELECTROSTRICTIVETENSORTOOLS_H

template <typename>
class RankFourTensorTempl;
typedef RankFourTensorTempl<Real> RankFourTensor;

namespace ElectrostrictiveTensorTools
{
  RankFourTensor computeProduct(const RankFourTensor & Cijkl, const RankFourTensor & Qmnkl);

  RankFourTensor computeProductQ(const RankFourTensor & Cijkl, const RankFourTensor & Qmnkl, const RankFourTensor & QQmnkl);

  /// Sum over (j,l) q_ijkl*v(j)*w(l)
  Real electrostrictiveProduct(const RankFourTensor & qijkl, unsigned int i, const RealVectorValue & v, unsigned int k, const RealVectorValue & w);

  /// Sum over l q_ijkl*v(j)
  Real electrostrictiveProduct(const RankFourTensor & qijkl, unsigned int i, const RealVectorValue & v,unsigned int k, const unsigned int l);
}

#endif //ELECTROSTRICTIVETENSORTOOLS_H
