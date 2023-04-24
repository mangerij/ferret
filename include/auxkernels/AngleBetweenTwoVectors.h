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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/
#ifndef ANGLEBETWEENTWOVECTORS_H
#define ANGLEBETWEENTWOVECTORS_H

#include "AuxKernel.h"

class AngleBetweenTwoVectors : public AuxKernel
{
public:
  AngleBetweenTwoVectors(const InputParameters & parameters);

  static InputParameters validParams();
  virtual ~AngleBetweenTwoVectors() {}

protected:
  virtual Real computeValue();

private:
  const VariableValue & _var1x;
  const VariableValue & _var1y;
  const VariableValue & _var1z;
  const VariableValue & _var2x;
  const VariableValue & _var2y;
  const VariableValue & _var2z;
  const Real _adjust_angle;

};

#endif
