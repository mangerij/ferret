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

#ifndef COMPUTEGCOEFFTENSOR_H
#define COMPUTEGCOEFFTENSOR_H

#include "RankFourTensor.h"
#include "ComputeRotatedGCoeffTensorBase.h"
#include "libmesh/quadrature.h"

class ComputeGCoeffTensor;

template<>
InputParameters validParams<ComputeGCoeffTensor>();

/**
 * ComputeGCoeffTensor defines an photostrictive tensor material object with a given base name.
 */
class ComputeGCoeffTensor : public ComputeRotatedGCoeffTensorBase
{
public:
  ComputeGCoeffTensor(const InputParameters & parameters);

protected:
  virtual void computeQpGCoeffTensor();

  /// Individual material information
  RankFourTensor _gijkl;
};

#endif //COMPUTEGCOEFFTENSOR_H
