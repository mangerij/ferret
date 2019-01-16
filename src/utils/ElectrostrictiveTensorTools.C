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

#include "MooseTypes.h"
#include "RankFourTensor.h"

#include "libmesh/vector_value.h"

namespace ElectrostrictiveTensorTools
{

RankFourTensor
computeProduct(const RankFourTensor & Cijkl, const RankFourTensor & Qmnkl)
{
  RankFourTensor result;
  ///Moose::out << "\n Performing C_ijmn Q_mnkl contraction on all the quadrature points?";
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j = 0; j < 3; ++j)
      for(unsigned int k = 0; k < 3; ++k)
        for(unsigned int l = 0; l < 3; ++l)
        {
          result(i,j,k,l) = 0.0;
          for(unsigned int m = 0; m < 3; ++m)
            for(unsigned int n = 0; n < 3; ++n)
            {
                ///sum += Cijkl(i, j, m, n) * Qmnkl(m, n, k, l);
                result(i,j,k,l) += Cijkl(i,j,m,n) * Qmnkl(m,n,k,l);
            }
            ///Moose::out << "\n q"; std::cout << i + 1 << j + 1 << k + 1 << l + 1; Moose::out << " = "; std::cout << result(i,j,k,l);
        }
  return result;
  ///Moose::out << "\n Complete.";
}

RankFourTensor

///here we use a different signature for the same member function
computeProductQ(const RankFourTensor & Cijkl, const RankFourTensor & Qmnkl, const RankFourTensor & QQmnkl)
{
  RankFourTensor result;
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j = 0; j < 3; ++j)
      for(unsigned int k = 0; k < 3; ++k)
        for(unsigned int l = 0; l < 3; ++l)
        {
          Real sum = 0.0;
          for(unsigned int m = 0; m < 3; ++m)
            for(unsigned int n = 0; n < 3; ++n)
              for(unsigned int r = 0; r < 3; ++r)
                for(unsigned int s = 0; s < 3; ++s)
                {
                  result(i,j,k,l) += Qmnkl(m, n, i, j) * Cijkl(m, n, r, s) * QQmnkl(r, s, k, l);
                }
        }
  return result;
}

Real
electrostrictiveProduct(const RankFourTensor & qijkl, unsigned int i, const RealVectorValue & v, unsigned int k, const RealVectorValue & w)
{
  // RankFourTensor qijkl;
  // Sum over (j,l) q_ijkl * v(j) * w(l) with k = _component
  Real sum = 0.0;
  for(unsigned int j = 0; j < 3; ++j)
    for(unsigned int l = 0; l < 3; ++l)
    {
      sum += qijkl(i, j, k, l) * v(j) * w(l);
    }
  return sum;
}

Real
electrostrictiveProduct(const RankFourTensor & qijkl, unsigned int i, const RealVectorValue & v, unsigned int k, const unsigned int l)
{
  // RankFourTensor qijkl;
  // Sum over j q_ijkl * v(j) where k and l = _component (used for DiagJacobian)
  Real sum = 0.0;
  for(unsigned int j = 0; j < 3; ++j)
    {
      sum += qijkl(i, j, k, l) * v(j);
    }
  return sum;
}

}
