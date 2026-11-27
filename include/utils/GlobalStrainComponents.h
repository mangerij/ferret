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

#pragma once

#include "MooseError.h"

#include <utility>

/**
 * Which strain component each entry of MOOSE's mean-strain SCALAR holds. The scalar's order is
 * tied to the mesh dimension by MOOSE's GlobalStrain kernel, and RankTwoTensor::fillFromScalarVariable
 * fills the strain as
 *   FIRST (1-D):  xx
 *   THIRD (2-D):  xx yy xy
 *   SIXTH (3-D):  xx yy zz yz xz xy
 * Kernels that work with the six 3-D components (indices 0..5 in the SIXTH order) look entries up here.
 */
namespace GlobalStrainComponents
{
/// entry b of a scalar of this order, as an index 0..5 into xx yy zz yz xz xy
inline unsigned int
sixIndex(unsigned int order, unsigned int b)
{
  static const unsigned int third[3] = {0, 1, 5};
  if (order == 1 && b == 0)
    return 0;
  if (order == 3 && b < 3)
    return third[b];
  if (order == 6 && b < 6)
    return b;
  mooseError("GlobalStrainComponents: the mean strain must be a FIRST, THIRD or SIXTH order SCALAR");
}

/// true for the shear entries (yz, xz, xy)
inline bool
isShear(unsigned int order, unsigned int b)
{
  return sixIndex(order, b) >= 3;
}

/// the tensor indices (k, l) of entry b
inline std::pair<unsigned int, unsigned int>
tensorIndices(unsigned int order, unsigned int b)
{
  static const std::pair<unsigned int, unsigned int> six[6] = {
      {0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}};
  return six[sixIndex(order, b)];
}
}
