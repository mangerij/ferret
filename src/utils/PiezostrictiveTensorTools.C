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
computeProduct(const RankFourTensor & Cpqij, const RankThreeTensor & Akpq)
{
  RankThreeTensor result;
  ///Moose::out << "\n Performing C_ijmn A_mnk contraction on all the quadrature points?";
  for (unsigned int k = 0; k < 3; ++k)
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        {
          result(k,i,j) = 0.0;
          for (unsigned int p = 0; p < 3; ++p)
            for (unsigned int q = 0; q < 3; ++q)
            {
                ///sum += Cijkl(i, j, m, n) * Amnk(m, n, k);
                result(k,i,j) += Cpqij(p,q,i,j) * Akpq(k,p,q);
                // Moose::out << "\n q"; std::cout << k + 1 << i + 1 << j + 1; Moose::out << " = "; std::cout << result(k,i,j);
            }
        }
  return result;
  ///Moose::out << "\n Complete.";
}


RankThreeTensor
computePiezoTransposeProduct(const RankFourTensor & Cijpq, const RankThreeTensor & Apqk)
{
  RankThreeTensor res;
  ///Moose::out << "\n Performing C_ijmn A_mnk contraction on all the quadrature points?";
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      for (unsigned int k = 0; k < 3; ++k)
        {
          res(i,j,k) = 0.0;
          for (unsigned int p = 0; p < 3; ++p)
            for (unsigned int q = 0; q < 3; ++q)
            {
                ///sum += Cijkl(i, j, m, n) * Amnk(m, n, k);
                res(i,j,k) += Cijpq(i,j,p,q) * Apqk(p,q,k);
                // Moose::out << "\n eT"; std::cout << i + 1 << j + 1 << k + 1; Moose::out << " = "; std::cout << res(i,j,k);
            }
        }
  return res;
  ///Moose::out << "\n Complete.";
}

Real
piezostrictiveProduct(const RankThreeTensor & Aijk, unsigned int i, const RealVectorValue & v, const RealVectorValue & w)
{
  /// RankThreeTensor Aijk;
  ///Sum over (j,l) A_ijk * v(j) * w(l) with k = _component
  Real sum = 0.0;
  for (unsigned int j = 0; j < 3; ++j)
    for (unsigned int k = 0; k < 3; ++k)
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
  for (unsigned int k = 0; k < 3; ++k)
    {
      sum += Aijk(i, j, k) * v(k);
    }
  return sum;
}

}
