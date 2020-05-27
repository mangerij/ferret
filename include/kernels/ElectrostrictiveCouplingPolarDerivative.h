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

#ifndef ELECTROSTRICTIVECOUPLINGDISTORTDERIVATIVE_H
#define ELECTROSTRICTIVECOUPLINGDISTORTDERIVATIVE_H

#include "Kernel.h"

//Forward Declarations
class ElectrostrictiveCouplingPolarDerivative;

template<>
InputParameters validParams<ElectrostrictiveCouplingPolarDerivative>();

class ElectrostrictiveCouplingPolarDerivative: public Kernel
{
public:
  ElectrostrictiveCouplingPolarDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _component;
  const unsigned int _disp_x_var;
  const unsigned int _disp_y_var;
  const unsigned int _disp_z_var;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableGradient & _disp_x_grad;
  const VariableGradient & _disp_y_grad;
  const VariableGradient & _disp_z_grad;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const MaterialProperty<Real> & _q11;
  const MaterialProperty<Real> & _q12;
  const MaterialProperty<Real> & _q44;
};
#endif //ELECTROSTRICTIVECOUPLINGDISTORTDERIVATIVE_H
