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

#ifndef DEMAGFIELDAUXPML_H
#define DEMAGFIELDAUXPML_H

#include "AuxKernel.h"

class DemagFieldAuxPML : public AuxKernel
{
public:
  DemagFieldAuxPML(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~DemagFieldAuxPML() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _component;
  const Real &   _deltasyminus;
  const Real &   _deltapyminus;
  const Real &   _deltawyminus;
  const Real &   _y0pmlminus;
  const VariableGradient & _phi1_grad;
  const VariableGradient & _potential_H_ext_grad;
};

#endif /* DEMAGFIELDAUXPML_H */
