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

#ifndef TIMEUSLL_H
#define TIMEUSLL_H

#include "TimeKernel.h"
#include "libmesh/quadrature.h"
#include "Assembly.h"

class TimeUSLL;

template<>
InputParameters validParams<TimeUSLL>();

class TimeUSLL : public TimeKernel
{
public:
  TimeUSLL(const InputParameters & parameters);

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

protected:

  const unsigned int _component;
  const unsigned int _azimuth_phi_var;
  const unsigned int _polar_theta_var;
  const VariableValue & _azimuth_phi;
  const VariableValue & _polar_theta;
  const VariableValue & _azimuth_phi_dot;
  const VariableValue & _polar_theta_dot;
  const VariableValue & _azimuth_phi_d_dot;
  const VariableValue & _polar_theta_d_dot;
};

#endif //TIMEUSLL_H
