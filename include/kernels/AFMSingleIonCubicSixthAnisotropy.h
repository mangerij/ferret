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

#ifndef AFMSINGLEIONCUBICSIXTHANISOTROPY_H
#define AFMSINGLEIONCUBICSIXTHANISOTROPY_H

#include "Kernel.h"

class AFMSingleIonCubicSixthAnisotropy: public Kernel
{
public:
  AFMSingleIonCubicSixthAnisotropy(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _mag_sub;
  const unsigned int _mag1_x_var;
  const unsigned int _mag1_y_var;
  const unsigned int _mag1_z_var;
  const unsigned int _mag2_x_var;
  const unsigned int _mag2_y_var;
  const unsigned int _mag2_z_var;
  const VariableValue & _mag1_x;
  const VariableValue & _mag1_y;
  const VariableValue & _mag1_z;
  const VariableValue & _mag2_x;
  const VariableValue & _mag2_y;
  const VariableValue & _mag2_z;
  const VariableValue & _antiphase_A_x;
  const VariableValue & _antiphase_A_y;
  const VariableValue & _antiphase_A_z;
  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> & _K1c;
  const MaterialProperty<Real> & _Kt;
  const MaterialProperty<Real> & _g0mu0Ms;
};
#endif //AFMSINGLEIONCUBICSIXTHANISOTROPY_H
