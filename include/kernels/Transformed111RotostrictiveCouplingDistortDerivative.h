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

#ifndef TRANSFORMED111ROTOSTRICTIVECOUPLINGDISTORTDERIVATIVE_H
#define TRANSFORMED111ROTOSTRICTIVECOUPLINGDISTORTDERIVATIVE_H

#include "Kernel.h"

class Transformed111RotostrictiveCouplingDistortDerivative: public Kernel
{
public:
  Transformed111RotostrictiveCouplingDistortDerivative(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _component;
  const unsigned int _u_x_var;
  const unsigned int _u_y_var;
  const unsigned int _u_z_var;
  const unsigned int _antiphase_A_x_var;
  const unsigned int _antiphase_A_y_var;
  const unsigned int _antiphase_A_z_var;
  const VariableGradient & _u_x_grad;
  const VariableGradient & _u_y_grad;
  const VariableGradient & _u_z_grad;
  const VariableValue & _antiphase_A_x;
  const VariableValue & _antiphase_A_y;
  const VariableValue & _antiphase_A_z;
  const MaterialProperty<Real> & _r11;
  const MaterialProperty<Real> & _r12;
  const MaterialProperty<Real> & _r44;
};
#endif //TRANSFORMED111ROTOSTRICTIVECOUPLINGDISTORTDERIVATIVE_H
