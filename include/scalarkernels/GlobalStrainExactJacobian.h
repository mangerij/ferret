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

#pragma once

#include "GlobalStrain.h"

/**
 * MOOSE's GlobalStrain scalar kernel with the exact eps_bar <- eps_bar Jacobian. MOOSE assembles
 * entry (b, c) as the single elasticity component C_(b)(c); a shear unit strain E_c moves both
 * eps_kl and eps_lk, so the true derivative of int sigma_(b) is C_(b)kl + C_(b)lk, twice MOOSE's
 * value on the three shear columns. Residual unchanged. Exact for translationally periodic cells
 * (as is MOOSE's own Jacobian: it ignores the non-periodic masking of ComputeGlobalStrain).
 */
class GlobalStrainExactJacobian : public GlobalStrain
{
public:
  static InputParameters validParams();

  GlobalStrainExactJacobian(const InputParameters & parameters);

  virtual void computeJacobian() override;
};
