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

#ifndef MAGHSTRONGCARTFECRCOALLOY_H
#define MAGHSTRONGCARTFECRCOALLOY_H

#include "Kernel.h"

class MagHStrongCartFeCrCoAlloy;

template<>
InputParameters validParams<MagHStrongCartFeCrCoAlloy>();

class MagHStrongCartFeCrCoAlloy: public Kernel
{
public:

  MagHStrongCartFeCrCoAlloy(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _c2_var;
  const unsigned int _c3_var;
  const VariableValue & _mag_x;
  const VariableValue & _mag_y;
  const VariableValue & _mag_z;
  const VariableValue & _c1;
  const VariableValue & _c2;
  const VariableValue & _c3;
  const Real _bohrM;
  const Real _T;

};
#endif //MAGHSTRONGCARTFECRCOALLOY_H