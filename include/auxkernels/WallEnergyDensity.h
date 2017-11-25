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

#ifndef WALLENERGYDENSITY_H
#define WALLENERGYDENSITY_H

#include "AuxKernel.h"


//Forward Declarations
class WallEnergyDensity;

template<>
InputParameters validParams<WallEnergyDensity>();

/**
 * Coupled auxiliary value
 */
class WallEnergyDensity : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  WallEnergyDensity(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  const VariableGradient & _polar_x_grad;
  const VariableGradient & _polar_y_grad;
  const VariableGradient & _polar_z_grad;
  const Real _G110,_G11, _G12, _G44, _G44P;
  const Real _len_scale;
};

#endif // WALLENERGYDENSITY_H
