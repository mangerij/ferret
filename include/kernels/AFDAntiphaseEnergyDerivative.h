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

#ifndef AFDANTIPHASEENERGYDERIVATIVE_H
#define AFDANTIPHASEENERGYDERIVATIVE_H

#include "Kernel.h"

class AFDAntiphaseEnergyDerivative;

template<>
InputParameters validParams<AFDAntiphaseEnergyDerivative>();

class AFDAntiphaseEnergyDerivative: public Kernel
{
public:

  AFDAntiphaseEnergyDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _antiferrodis_A_x_var;
  const unsigned int _antiferrodis_A_y_var;
  const unsigned int _antiferrodis_A_z_var;
  const VariableGradient & _antiferrodis_A_i_grad;
  const VariableGradient & _antiferrodis_A_j_grad;
  const VariableGradient & _antiferrodis_A_k_grad;
  const unsigned int _ii, _jj, _kk;
  const Real _H110, _H11, _H12, _H44, _H44P;
  const Real _len_scale;

};
#endif //AFDANTIPHASEENERGYDERIVATIVE_H
