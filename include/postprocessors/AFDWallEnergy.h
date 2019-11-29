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

#ifndef AFDWALLENERGY_H
#define AFDWALLENERGY_H

#include "ElementIntegralPostprocessor.h"

class AFDWallEnergy;

template<>
InputParameters validParams<AFDWallEnergy>();

class AFDWallEnergy : public ElementIntegralPostprocessor
{
public:
  AFDWallEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableGradient& _antiferrodis_A_x_grad;
  const VariableGradient& _antiferrodis_A_y_grad;
  const VariableGradient& _antiferrodis_A_z_grad;
  const Real _H110,_H11, _H12, _H44, _H44P;
  const Real _len_scale;
};

#endif
