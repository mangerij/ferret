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

#ifndef INDUCEDPWALLENERGYDENSITY_H
#define INDUCEDPWALLENERGYDENSITY_H

#include "AuxKernel.h"

class InducedPWallEnergyDensity;

template<>
InputParameters validParams<InducedPWallEnergyDensity>();

class InducedPWallEnergyDensity : public AuxKernel
{
public:
  InducedPWallEnergyDensity(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  const VariableGradient & _polar_x_grad;
  const VariableGradient & _polar_y_grad;
  const VariableGradient & _polar_z_grad;
  const VariableGradient & _induced_polar_x_grad;
  const VariableGradient & _induced_polar_y_grad;
  const VariableGradient & _induced_polar_z_grad;
  const Real _G110,_G11, _G12, _G44, _G44P;
  const Real _len_scale;
};

#endif // INDUCEDPWALLENERGYDENSITY_H
