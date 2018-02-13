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
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
//#include "RankSixTensor.h"
//#include "RankEightTensor.h"

namespace LandauTensorTools
{

//the below are useful for Postprocessors and AuxKernels that define the various free energy densities.

Real
landauTwoProduct(const RankTwoTensor & aij, const RealVectorValue & p)
{
  Real sum = 0.0;
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j = 0; j < 3; ++j)
    {
      sum += aij(i, j) * p(i) * p(j);
    }
  return sum;
}

Real
landauFourProduct(const RankFourTensor & aijkl, const RealVectorValue & p)
{
  Real sum = 0.0;
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j = 0; j < 3; ++j)
      for(unsigned int k = 0; k < 3; ++k)
        for(unsigned int l = 0; l < 3; ++l)
        {
          sum += aijkl(i, j, k, l) * p(i) * p(j) * p(k) * p(l);
        }
  return sum;
}

//Real
//landauSixProduct(const RankSixTensor & aijklmn, const RealVectorValue & p)
//{
//  Real sum = 0.0;
//  for(unsigned int i = 0; i < 3; ++i)
//    for(unsigned int j = 0; j < 3; ++j)
//      for(unsigned int k = 0; k < 3; ++k)
//        for(unsigned int l = 0; l < 3; ++l)
//          for(unsigned int m = 0; m < 3; ++m)
//            for(unsigned int n = 0; n < 3; ++n)
//            {
//              sum += aijklmn(i, j, k, l, m, n) * p(i) * p(j) * p(k) * p(l) * p(m) * p(n);
//            }
//  return sum;
//}

//Real
//landauEightProduct(const RankEightTensor & aijklmnrs, const RealVectorValue & p)
//{
//  Real sum = 0.0;
//  for(unsigned int i = 0; i < 3; ++i)
//    for(unsigned int j = 0; j < 3; ++j)
//      for(unsigned int k = 0; k < 3; ++k)
//        for(unsigned int l = 0; l < 3; ++l)
//          for(unsigned int m = 0; m < 3; ++m)
//            for(unsigned int n = 0; n < 3; ++n)
//              for(unsigned int r = 0; r < 3; ++r)
//                for(unsigned int s = 0; s < 3; ++s)
//                {
//                  sum += aijklmnrs(i, j, k, l, m, n, r, s) * p(i) * p(j) * p(k) * p(l) * p(m) * p(n) * p(r) * p(s);
//                }
//  return sum;
//}

//the below are useful for Kernels that define the variations of the free energy contributions.

Real
landauTwoProductDerivative(const RankTwoTensor & aij, unsigned int k, const RealVectorValue & p)
{
  Real sum = 0.0;
  for(unsigned int j = 0; j < 3; ++j)
  {
    sum += aij(k, j) * p(j) + aij(j, k) * p(j);
  }
  return sum;
}

Real
landauFourProductDerivative(const RankFourTensor & aijkl, unsigned int m, const RealVectorValue & p)
{
  Real sum = 0.0;
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j = 0; j < 3; ++j)
      for(unsigned int k = 0; k < 3; ++k)
      {
        sum += aijkl(i, j, k, m) * p(i) * p(j) * p(k) + aijkl(i, j, m, k) * p(i) * p(j) * p(k) + aijkl(i, m, j, k) * p(i) * p(j) * p(k) + aijkl(m, i, j, k) * p(i) * p(j) * p(k);
      }
  return sum;
}

//need jacobians and off diagonal jacobians :(


}
