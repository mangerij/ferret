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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#ifndef LOCALSATURATIONCONSTRAINT_H
#define LOCALSATURATIONCONSTRAINT_H

#include "Kernel.h"

class LocalSaturationConstraint;

template<>
InputParameters validParams<LocalSaturationConstraint>();

class  LocalSaturationConstraint: public Kernel
{
public:
  LocalSaturationConstraint(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  const unsigned int _component;
  const unsigned int _mag_x_var;
  const unsigned int _mag_y_var;
  const unsigned int _mag_z_var;
  const VariableValue & _mag_x;
  const VariableValue & _mag_y;
  const VariableValue & _mag_z;
  const Real _var_mag;
  const Real _sat_penalty;

};
#endif //LOCALSATURATIONCONSTRAINT_H
