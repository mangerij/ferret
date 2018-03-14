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

#ifndef CONSTRAINEDINTERACTIONLLG_H
#define CONSTRAINEDINTERACTIONLLG_H

#include "Kernel.h"

class ConstrainedInteractionLLG;

template<>
InputParameters validParams<ConstrainedInteractionLLG>();

class ConstrainedInteractionLLG: public Kernel
{
public:

  ConstrainedInteractionLLG(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _component;
  const unsigned int _potential_H_int_var;
  const unsigned int _potential_H_ext_var;
  const VariableGradient & _potential_H_int_grad;
  const VariableGradient & _potential_H_ext_grad;
  const unsigned int _azimuth_phi_var;
  const unsigned int _polar_theta_var;
  const VariableValue & _azimuth_phi;
  const VariableValue & _polar_theta;
  const Real _alpha;
  const Real _mu0;
  const Real _g0;
  const Real _M;
};
#endif //CONSTRAINEDINTERACTIONLLG_H
