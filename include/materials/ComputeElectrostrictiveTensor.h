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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#ifndef COMPUTEELECTROSTRICTIVETENSOR_H
#define COMPUTEELECTROSTRICTIVETENSOR_H

#include "RankFourTensor.h"
#include "ElectrostrictiveTensorTools.h"
#include "ComputeRotatedElectrostrictiveTensorBase.h"
#include "libmesh/quadrature.h"

class ComputeElectrostrictiveTensor;

template<>
InputParameters validParams<ComputeElectrostrictiveTensor>();

/**
 * ComputeElectrostrictiveTensor defines an electrostrictive tensor material object with a given base name.
 */
class ComputeElectrostrictiveTensor : public ComputeRotatedElectrostrictiveTensorBase
{
public:
  ComputeElectrostrictiveTensor(const InputParameters & parameters);

protected:
  virtual void computeQpElectrostrictiveTensor();

  /// Individual material information
  bool _compute_electrostrictive_coeff;
  RankFourTensor _Qmnkl;
  RankFourTensor _Cijkl;
  RankFourTensor _qijkl;

 private:
   const MaterialProperty<RankFourTensor> & _elasticity_tensor;
};

#endif //COMPUTEELECTROSTRICTIVETENSOR_H
