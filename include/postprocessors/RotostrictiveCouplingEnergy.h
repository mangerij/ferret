/**
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
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#ifndef ROTOSTRICTIVECOUPLINGENERGY_H
#define ROTOSTRICTIVECOUPLINGENERGY_H

#include "ElementIntegralPostprocessor.h"

class RotostrictiveCouplingEnergy;

template<>
InputParameters validParams<RotostrictiveCouplingEnergy>();


class RotostrictiveCouplingEnergy : public ElementIntegralPostprocessor
{
public:
  RotostrictiveCouplingEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

private:
  const unsigned int _disp_x_var;
  const unsigned int _disp_y_var;
  const unsigned int _disp_z_var;
  const unsigned int _antiferrodis_A_x_var;
  const unsigned int _antiferrodis_A_y_var;
  const unsigned int _antiferrodis_A_z_var;
  const VariableGradient & _disp_x_grad;
  const VariableGradient & _disp_y_grad;
  const VariableGradient & _disp_z_grad;
  const VariableValue & _antiferrodis_A_x;
  const VariableValue & _antiferrodis_A_y;
  const VariableValue & _antiferrodis_A_z;
  const Real _r11;
  const Real _r12;
  const Real _r44;
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm
};

#endif
