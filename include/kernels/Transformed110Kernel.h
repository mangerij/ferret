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

#ifndef TRANSFORMED110KERNEL_H
#define TRANSFORMED110KERNEL_H

#include "Kernel.h"

class Transformed110Kernel: public Kernel
{
public:

  Transformed110Kernel(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _order_param_x_var;
  const unsigned int _order_param_y_var;
  const unsigned int _order_param_z_var;
  const VariableValue & _fb_x;
  const VariableValue & _fb_y;
  const VariableValue & _fb_z;
  const VariableValue & _Jb_xx;
  const VariableValue & _Jb_yy;
  const VariableValue & _Jb_zz;
  const VariableValue & _Jb_xy;
  const VariableValue & _Jb_yz;
  const VariableValue & _Jb_xz;

};
#endif //TRANSFORMED110KERNEL_H
