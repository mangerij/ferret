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

#ifndef EXCHANGECARTLL_H
#define EXCHANGECARTLL_H

#include "Kernel.h"

class  ExchangeCartLL: public Kernel
{
public:
  ExchangeCartLL(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _mag_x_var;
  const unsigned int _mag_y_var;
  const unsigned int _mag_z_var;
  const VariableValue & _mag_x;
  const VariableValue & _mag_y;
  const VariableValue & _mag_z;
  const VariableGradient & _mag_x_grad;
  const VariableGradient & _mag_y_grad;
  const VariableGradient & _mag_z_grad;
  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> & _g0;
  const MaterialProperty<Real> & _Ae;
  const MaterialProperty<Real> & _Ms;

};
#endif //EXCHANGECARTLL_H
