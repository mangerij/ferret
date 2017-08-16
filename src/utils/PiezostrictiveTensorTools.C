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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "MooseTypes.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"

namespace PiezostrictiveTensorTools
{

RankThreeTensor
computeProduct(const RankFourTensor & Cijmn, const RankThreeTensor & Amnk)
{
  RankThreeTensor result;
  ///Moose::out << "\n Performing C_ijmn A_mnk contraction on all the quadrature points?";
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j = 0; j < 3; ++j)
      for(unsigned int k = 0; k < 3; ++k)
        {
          result(i,j,k) = 0.0;
          for(unsigned int m = 0; m < 3; ++m)
            for(unsigned int n = 0; n < 3; ++n)
            {
                ///sum += Cijkl(i, j, m, n) * Amnk(m, n, k);
                result(i,j,k) += Cijmn(i,j,m,n) * Amnk(m,n,k);
            }
            ///Moose::out << "\n q"; std::cout << i + 1 << j + 1 << k + 1 << l + 1; Moose::out << " = "; std::cout << result(i,j,k,l);
        }
  return result;
  ///Moose::out << "\n Complete.";
}


Real
piezostrictiveProduct(const RankThreeTensor & Aijk, unsigned int i, const RealVectorValue & v, const RealVectorValue & w)
{
  /// RankThreeTensor Aijk;
  ///Sum over (j,l) A_ijk * v(j) * w(l) with k = _component
  Real sum = 0.0;
  for(unsigned int j = 0; j < 3; ++j)
    for(unsigned int k = 0; k < 3; ++k)
    {
      sum += Aijk(i, j, k) * v(j) * w(k);
    }
  return sum;
}

Real
piezostrictiveProduct(const RankThreeTensor & Aijk, unsigned int i, const RealVectorValue & v, unsigned int j)
{
  /// RankFourTensor qijkl;
  ///Sum over j q_ijkl * v(j) where k and l = _component (used for DiagJacobian)
  Real sum = 0.0;
  for(unsigned int k = 0; k < 3; ++k)
    {
      sum += Aijk(i, j, k) * v(k);
    }
  return sum;
}

}
