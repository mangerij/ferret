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

#ifndef COMPUTEROTATEDINDICATRIXBASE_H
#define COMPUTEROTATEDINDICATRIXBASE_H

#include "ComputeIndicatrixBase.h"
#include "RankTwoTensor.h"

/**
 * ComputeRotatedIndicatrixBase the base class for computing photostrictive tensors
 */
class ComputeRotatedIndicatrixBase : public ComputeIndicatrixBase
{
public:
  ComputeRotatedIndicatrixBase(const InputParameters & parameters);

protected:
  RealVectorValue _Euler_angles;
};

#endif //COMPUTEROTATEDINDICATRIXBASE_H
