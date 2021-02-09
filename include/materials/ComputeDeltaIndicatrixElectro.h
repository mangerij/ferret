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

#ifndef COMPUTEDELTAINDICATRIXELECTRO_H
#define COMPUTEDELTAINDICATRIXELECTRO_H

#include "RankThreeTensor.h"
#include "ComputeDeltaIndicatrixElectroBase.h"
#include "libmesh/quadrature.h"

class ComputeDeltaIndicatrixElectro;

template<>
InputParameters validParams<ComputeDeltaIndicatrixElectro>();

/**
 * ComputeDeltaIndicatrixElectro defines an impermeability tensor material object with a given base name.
 */
class ComputeDeltaIndicatrixElectro : public ComputeDeltaIndicatrixElectroBase
{
public:
  ComputeDeltaIndicatrixElectro(const InputParameters & parameters);

protected:
  virtual void computeQpDeltaIndicatrixElectro();
  const MaterialProperty<RankThreeTensor> & _electrooptic_tensor;
  const VariableGradient & _potential_E_int_grad;

};

#endif //COMPUTEDELTAINDICATRIXELECTRO_H
