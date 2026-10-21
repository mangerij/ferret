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

#ifndef MASTERZLSPINTRANSFERTORQUECARTLLG_H
#define MASTERZLSPINTRANSFERTORQUECARTLLG_H

#include "Kernel.h"

class MasterZLSpinTransferTorqueCartLLG : public Kernel
{
public:
  static InputParameters validParams();

  MasterZLSpinTransferTorqueCartLLG(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  Real torqueComponent(unsigned int component) const;
  Real torqueDerivative(unsigned int component, unsigned int variable_component) const;
  int leviCivita(unsigned int i, unsigned int j, unsigned int k) const;

  /// not const: TangentPlaneProjected evaluates all three components through it
  unsigned int _component;
  const unsigned int _current_direction;

  const unsigned int _mag_x_var;
  const unsigned int _mag_y_var;
  const unsigned int _mag_z_var;

  const VariableValue & _mag_x;
  const VariableValue & _mag_y;
  const VariableValue & _mag_z;

  const VariableGradient & _mag_x_grad;
  const VariableGradient & _mag_y_grad;
  const VariableGradient & _mag_z_grad;

  const MaterialProperty<Real> & _alpha;

  const Real _u;
  const Real _xi;
};

#endif
