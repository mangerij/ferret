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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#pragma once

#include "ODEKernel.h"

class ScalarElasticEnergy : public ODEKernel
{
public:

  ScalarElasticEnergy(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobianScalar(unsigned int jvar) override;

  const unsigned int _component;

  unsigned int _e_xx_var;
  unsigned int _e_yy_var;
  unsigned int _e_zz_var;
  unsigned int _e_xy_var;
  unsigned int _e_yz_var;
  unsigned int _e_zx_var;

  const VariableValue & _e_xx;
  const VariableValue & _e_yy;
  const VariableValue & _e_zz;
  const VariableValue & _e_xy;
  const VariableValue & _e_yz;
  const VariableValue & _e_zx;

  const Real _G;
  const Real _C11;
  const Real _C12;
  const Real _C44;

};
