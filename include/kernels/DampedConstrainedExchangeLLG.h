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

#ifndef DAMPEDCONSTRAINEDEXCHANGELLG_H
#define DAMPEDCONSTRAINEDEXCHANGELLG_H

#include "Kernel.h"

class DampedConstrainedExchangeLLG;

template<>
InputParameters validParams<DampedConstrainedExchangeLLG>();

class  DampedConstrainedExchangeLLG: public Kernel
{
public:
  DampedConstrainedExchangeLLG(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _azimuth_phi_var;
  const unsigned int _polar_theta_var;
  const VariableValue & _azimuth_phi;
  const VariableValue & _polar_theta;
  const VariableGradient & _azimuth_phi_grad;
  const VariableGradient & _polar_theta_grad;
  const Real _alpha;
  const Real _g0;
  const Real _Ae;
  const Real _M;


};
#endif //DAMPEDCONSTRAINEDEXCHANGELLG_H
