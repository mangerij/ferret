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

#ifndef COMPUTEROTATEDRANKSIXLANDAUTENSORBASE_H
#define COMPUTEROTATEDRANKSIXLANDAUTENSORBASE_H

#include "ComputeRankSixLandauTensorBase.h"

class ComputeRotatedRankSixLandauTensorBase;

template<>
InputParameters validParams<ComputeRotatedRankSixLandauTensorBase>();

/**
 * ComputeRotatedRankSixLandauTensorBase is an intermediate base class that rotates the \alpha_{ijklmn} tensor based on euler angles.
 */
class ComputeRotatedRankSixLandauTensorBase : public ComputeRankSixLandauTensorBase
{
public:
  ComputeRotatedRankSixLandauTensorBase(const InputParameters & parameters);

protected:
  RealVectorValue _Euler_angles;
};

#endif //COMPUTEROTATEDRANKSIXLANDAUTENSORBASE_H