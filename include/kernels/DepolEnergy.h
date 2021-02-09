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

#ifndef DEPOLENERGY_H
#define DEPOLENERGY_H

#include "Kernel.h"

class DepolEnergy;

template<>
InputParameters validParams<DepolEnergy>();

class DepolEnergy: public Kernel
{
public:

  DepolEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const unsigned int _polar_z_var;
  const VariableValue & _polar_z;

  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm
  const PostprocessorValue & _avePz;
  const Real _lambda;
  const Real _permitivitty;

};
#endif //DEPOLENERGY_H
