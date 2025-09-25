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

#ifndef QUASISTATICFIELDAUX_H
#define QUASISTATICFIELDAUX_H

#include "AuxKernel.h"

class QuasistaticFieldAux : public AuxKernel
{
public:
  QuasistaticFieldAux(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~QuasistaticFieldAux() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _component;
  const VariableGradient & _potential_int_grad;
  const VariableGradient & _potential_ext_grad;
  const Real _conversion_factor;
};

#endif /* QUASISTATICFIELDAUX_H */
