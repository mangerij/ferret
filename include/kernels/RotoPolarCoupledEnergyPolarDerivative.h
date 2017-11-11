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

#ifndef ROTOPOLARCOUPLEDENERGYPOLARDERIVATIVE_H
#define ROTOPOLARCOUPLEDENERGYPOLARDERIVATIVE_H

#include "Kernel.h"

class RotoPolarCoupledEnergyPolarDerivative;

template<>
InputParameters validParams<RotoPolarCoupledEnergyPolarDerivative>();

class RotoPolarCoupledEnergyPolarDerivative: public Kernel
{
public:

  RotoPolarCoupledEnergyPolarDerivative(const InputParameters & parameters);

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
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const Real _t1111, _t1122, _t1212,  _t42111111, _t24111111, _t42111122,  _t24112222, _t42112233, _t24112233, _t42112211, _t24111122, _t42111212, _t42123312, _t24121112, _t24121233, _t6211111111, _t2611111111, _t6211111122, _t2611222222, _t4411111111, _t4411112222;
  const Real _len_scale;
};
#endif //ROTOPOLARCOUPLEDENERGYPOLARDERIVATIVE_H
