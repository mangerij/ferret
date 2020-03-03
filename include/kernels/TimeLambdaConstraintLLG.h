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

#ifndef TIMELAMBDACONSTRAINTLLG_H
#define TIMELAMBDACONSTRAINTLLG_H

#include "TimeKernel.h"
#include "libmesh/quadrature.h"
#include "Assembly.h"

class TimeLambdaConstraintLLG;

template<>
InputParameters validParams<TimeLambdaConstraintLLG>();

class TimeLambdaConstraintLLG : public TimeKernel
{
public:
  TimeLambdaConstraintLLG(const InputParameters & parameters);

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

protected:

  const unsigned int _mag_x_var;
  const unsigned int _mag_y_var;
  const unsigned int _mag_z_var;
  const unsigned int _lambda_var;

  const VariableValue & _mag_x;
  const VariableValue & _mag_y;
  const VariableValue & _mag_z;
  const VariableValue & _lambda;

  const VariableValue & _mag_x_dot;
  const VariableValue & _mag_y_dot;
  const VariableValue & _mag_z_dot;

  const VariableValue & _mag_x_d_dot;
  const VariableValue & _mag_y_d_dot;
  const VariableValue & _mag_z_d_dot;

  const Real _eps;
};

#endif //TIMELAMBDACONSTRAINTLLG_H
