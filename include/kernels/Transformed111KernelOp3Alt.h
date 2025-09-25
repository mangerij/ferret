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

#ifndef TRANSFORMED111KERNELOP3ALT_H
#define TRANSFORMED111KERNELOP3ALT_H

#include "Kernel.h"

class Transformed111KernelOp3Alt: public Kernel
{
public:

  Transformed111KernelOp3Alt(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _order_param_x_var;
  const unsigned int _order_param_y_var;
  const unsigned int _order_param_z_var;
  const unsigned int _order_param2_x_var;
  const unsigned int _order_param2_y_var;
  const unsigned int _order_param2_z_var;

  const VariableValue & _f_q0;
  const VariableValue & _f_q1;
  const VariableValue & _f_q2;

  const VariableValue & _J_q0q0;
  const VariableValue & _J_q1q1;
  const VariableValue & _J_q2q2;

  const VariableValue & _J_q0q1;
  const VariableValue & _J_q0q2;
  const VariableValue & _J_q0q3;
  const VariableValue & _J_q0q4;
  const VariableValue & _J_q0q5;

  const VariableValue & _J_q1q0;
  const VariableValue & _J_q1q2;
  const VariableValue & _J_q1q3;
  const VariableValue & _J_q1q4;
  const VariableValue & _J_q1q5;

  const VariableValue & _J_q2q0;
  const VariableValue & _J_q2q1;
  const VariableValue & _J_q2q3;
  const VariableValue & _J_q2q4;
  const VariableValue & _J_q2q5;

};
#endif //TRANSFORMED111KERNELOP3ALT_H
