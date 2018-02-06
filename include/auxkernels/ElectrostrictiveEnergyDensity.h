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
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#ifndef ELECTROSTRICTIVEENERGYDENSITY_H
#define ELECTROSTRICTIVEENERGYDENSITY_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "ComputeElectrostrictiveTensor.h"
#include "ComputeEigenstrain.h"

//Forward Declarations
class ElectrostrictiveEnergyDensity;

template<>
InputParameters validParams<ElectrostrictiveEnergyDensity>();

/**
 * Coupled auxiliary value
 */
class ElectrostrictiveEnergyDensity : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  ElectrostrictiveEnergyDensity(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  const MaterialProperty<RankFourTensor> & _electrostrictive_tensor;
  const MaterialProperty<RankTwoTensor> & _eigenstrain;
  const VariableGradient & _disp_x_grad;
  const VariableGradient & _disp_y_grad;
  const VariableGradient & _disp_z_grad;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const Real _artificial;
  const Real _len_scale;
};

#endif // ELECTROSTRICTIVEENERGYDENSITY_H
