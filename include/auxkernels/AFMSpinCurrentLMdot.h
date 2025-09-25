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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#ifndef AFMSPINCURRENTLMDOT_H
#define AFMSPINCURRENTLMDOT_H

#include "AuxKernel.h"

class AFMSpinCurrentLMdot : public AuxKernel
{
public:
  AFMSpinCurrentLMdot(const InputParameters & parameters);

  static InputParameters validParams();
  virtual ~AFMSpinCurrentLMdot() {}

protected:
  virtual Real computeValue();

private:
  const VariableValue & _Neel_L_x;
  const VariableValue & _Neel_L_y;
  const VariableValue & _Neel_L_z;
  const VariableValue & _dSSmag_dt_x;
  const VariableValue & _dSSmag_dt_y;
  const VariableValue & _dSSmag_dt_z;
  const unsigned int _component;
  const Real _factor;
};

#endif

