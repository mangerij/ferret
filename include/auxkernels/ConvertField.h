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

#ifndef CONVERTFIELD_H
#define CONVERTFIELD_H

#include "AuxKernel.h"

class ConvertField : public AuxKernel
{
public:
  ConvertField(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~ConvertField() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _conv_type;
  const unsigned int _component;
  const VariableValue & _polar_x_red;
  const VariableValue & _polar_y_red;
  const VariableValue & _polar_z_red;
  const VariableValue & _antiferrodis_A_x_red;
  const VariableValue & _antiferrodis_A_y_red;
  const VariableValue & _antiferrodis_A_z_red;
  const VariableValue & _disp_x_red;
  const VariableValue & _disp_y_red;
  const VariableValue & _disp_z_red;
  const VariableValue & _strain_xx_red;
  const VariableValue & _strain_yy_red;
  const VariableValue & _strain_zz_red;
  const VariableValue & _strain_xy_red;
  const VariableValue & _strain_xz_red;
  const VariableValue & _strain_yz_red;
  const Real _Ps, _As, _uscale, _es_norm, _es_shear;
};

#endif // CONVERTFIELD_H
