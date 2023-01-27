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

#ifndef BULKENERGYDERIVATIVESIXTHCOUPLEDT_H
#define BULKENERGYDERIVATIVESIXTHCOUPLEDT_H

#include "Kernel.h"

class BulkEnergyDerivativeSixthCoupledT: public Kernel
{
public:

  BulkEnergyDerivativeSixthCoupledT(const InputParameters & parameters);

  static InputParameters validParams();

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
  const unsigned int _temperature_var;
  const VariableValue & _temperature;
  const Real _alpha0, _alpha11, _alpha12, _alpha111, _alpha112, _alpha123, _Tc;
  const Real _len_scale;
};
#endif //BULKENERGYDERIVATIVESIXTHCOUPLEDT_H
