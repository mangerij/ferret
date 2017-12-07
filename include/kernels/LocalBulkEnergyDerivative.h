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

#ifndef LOCALBULKENERGYDERIVATIVE_H
#define LOCALBULKENERGYDERIVATIVE_H

#include "Kernel.h"

class LocalBulkEnergyDerivative;

template<>
InputParameters validParams<LocalBulkEnergyDerivative>();

class LocalBulkEnergyDerivative: public Kernel
{
public:

  LocalBulkEnergyDerivative(const InputParameters & parameters);

  static constexpr Real _default_val = 123456.0;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  bool _isRenorm;
  const Real _alpha1, _alpha3, _alpha11, _alpha33, _alpha12, _alpha13, _alpha111, _alpha112, _alpha123, _alpha1111, _alpha1112, _alpha1122, _alpha1123;
  const Real _len_scale;
};
#endif //LOCALBULKENERGYDERIVATIVE_H
