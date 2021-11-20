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

class ScalarRotopolarEnergy : public ODEKernel
{
public:
  /**
   * Constructor
   */
  ScalarRotopolarEnergy(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobianScalar(unsigned int jvar) override;

  const unsigned int _component;

  unsigned int _polar_x_var;
  unsigned int _polar_y_var;
  unsigned int _polar_z_var;

  unsigned int _antiferrodis_A_x_var;
  unsigned int _antiferrodis_A_y_var;
  unsigned int _antiferrodis_A_z_var;

  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;

  const VariableValue & _antiferrodis_A_x;
  const VariableValue & _antiferrodis_A_y;
  const VariableValue & _antiferrodis_A_z;

  const Real _G;
  const Real _Bpr;
  const Real _Cpr;
  const Real _Cprp;
};
