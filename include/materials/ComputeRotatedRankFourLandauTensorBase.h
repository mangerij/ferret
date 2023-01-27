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

#ifndef COMPUTEROTATEDRANKFOURLANDAUTENSORBASE_H
#define COMPUTEROTATEDRANKFOURLANDAUTENSORBASE_H

#include "ComputeRankFourLandauTensorBase.h"

/**
 * ComputeRotatedRankFourLandauTensorBase is an intermediate base class that rotates the \alpha_{ij} tensor based on euler angles.
 */
class ComputeRotatedRankFourLandauTensorBase : public ComputeRankFourLandauTensorBase
{
public:
  ComputeRotatedRankFourLandauTensorBase(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  RealVectorValue _Euler_angles;
};

#endif //COMPUTEROTATEDRANKFOURLANDAUTENSORBASE_H
