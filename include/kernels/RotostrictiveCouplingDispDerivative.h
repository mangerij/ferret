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

#ifndef ROTOSTRICTIVECOUPLINGDISPDERIVATIVE_H
#define ROTOSTRICTIVECOUPLINGDISPDERIVATIVE_H

#include "Kernel.h"

//Forward Declarations
class RotostrictiveCouplingDispDerivative;

template<>
InputParameters validParams<RotostrictiveCouplingDispDerivative>();

class RotostrictiveCouplingDispDerivative: public Kernel
{
public:

  RotostrictiveCouplingDispDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _component;
  const unsigned int _antiferrodis_A_x_var;
  const unsigned int _antiferrodis_A_y_var;
  const unsigned int _antiferrodis_A_z_var;
  const VariableValue & _antiferrodis_A_x;
  const VariableValue & _antiferrodis_A_y;
  const VariableValue & _antiferrodis_A_z;
  const Real _C11;
  const Real _C12;
  const Real _C44;
  const Real _R11;
  const Real _R12;
  const Real _R44;
  const Real _len_scale; //dimension unit, eg: 1e-9 for nm

};
#endif //ROTOSTRICTIVECOUPLINGDISPDERIVATIVE_H