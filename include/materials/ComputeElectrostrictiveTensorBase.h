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

#ifndef COMPUTEELECTROSTRICTIVETENSORBASE_H
#define COMPUTEELECTROSTRICTIVETENSORBASE_H

#include "Material.h"
#include "RankFourTensor.h"

class ComputeElectrostrictiveTensorBase;

template<>
InputParameters validParams<ComputeElectrostrictiveTensorBase>();

/**
 * ComputeElectrostrictiveTensorBase the base class for computing electrostrictive tensors
 */
class ComputeElectrostrictiveTensorBase : public Material
{
public:
  ComputeElectrostrictiveTensorBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpElectrostrictiveTensor() = 0;

  std::string _base_name;
  std::string _electrostrictive_tensor_name;
  std::string _electrostrictive_coefficients_name;

  MaterialProperty<RankFourTensor> & _electrostrictive_tensor;
  MaterialProperty<RankFourTensor> & _electrostrictive_coefficients;
};

#endif //COMPUTEELECTROSTRICTIVETENSORBASE_H
