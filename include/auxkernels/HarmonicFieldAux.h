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

#ifndef HARMONICFIELDAUX_H
#define HARMONICFIELDAUX_H

#include "AuxKernel.h"

class HarmonicFieldAux : public AuxKernel
{
public:
  HarmonicFieldAux(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~HarmonicFieldAux() {}

protected:
  virtual Real computeValue();

private:
  const Real _amplitude;
  const Real _correction;
  const Real _frequency;
  const Real _tshift;
  const Real _ton;
};

#endif // HARMONICFIELDAUX_H