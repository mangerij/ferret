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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#ifndef HEATFLOWELECTRICT_H
#define HEATFLOWELECTRICT_H

#include "Kernel.h"

class HeatFlowElectricT : public Kernel
{
public:
  HeatFlowElectricT(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _potential_E_int_var;
  const VariableValue & _potential_E_int;
  const VariableGradient & _potential_E_int_grad;
  const unsigned int _T_var;
  const VariableValue & _T;
  const VariableGradient & _T_grad;
  // const Real _electrical_conductivity;
  // const Real _seebeck_coefficient;
  // const Real _thermal_conductivity;
  const MaterialProperty<Real> & _ecC;
  const MaterialProperty<Real> & _sbC;
  const MaterialProperty<Real> & _thC;
};
#endif
