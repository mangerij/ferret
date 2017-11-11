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

#ifndef ROTOBULKENERGYDERIVATIVEEIGHTH_H
#define ROTOBULKENERGYDERIVATIVEEIGHTH_H

#include "Kernel.h"

class RotoBulkEnergyDerivativeEighth;

template<>
InputParameters validParams<RotoBulkEnergyDerivativeEighth>();

class RotoBulkEnergyDerivativeEighth: public Kernel
{
public:

  RotoBulkEnergyDerivativeEighth(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _antiferrodis_A_x_var;
  const unsigned int _antiferrodis_A_y_var;
  const unsigned int _antiferrodis_A_z_var;
  const VariableValue & _antiferrodis_A_x;
  const VariableValue & _antiferrodis_A_y;
  const VariableValue & _antiferrodis_A_z;
  const Real _beta1, _beta11, _beta12, _beta111, _beta112,_beta123, _beta1111, _beta1112, _beta1122, _beta1123;
  const Real _len_scale;
};
#endif //ROTOBULKENERGYDERIVATIVEEIGHTH_H
