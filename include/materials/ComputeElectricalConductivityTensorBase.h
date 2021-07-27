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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#ifndef COMPUTEELECTRICALCONDUCTIVITYTENSORBASE_H
#define COMPUTEELECTRICALCONDUCTIVITYTENSORBASE_H

#include "Material.h"
#include "RankTwoTensor.h"

class ComputeElectricalConductivityTensorBase;

template <>
InputParameters validParams<ComputeElectricalConductivityTensorBase>();

/**
 * ComputeElectricalConductivityTensorBase the base class for computing thermal conductivity tensors
 */
class ComputeElectricalConductivityTensorBase : public Material
{
public:
  ComputeElectricalConductivityTensorBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpElectricalConductivityTensor() = 0;

  std::string _base_name;
  std::string _ecC_tensor_name;

  MaterialProperty<RankTwoTensor> & _ecC_tensor;
};

#endif // COMPUTEElectricalConductivityTENSORBASE_H
