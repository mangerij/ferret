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

#ifndef GLOBALSTRAINOFFDIAG_H
#define GLOBALSTRAINOFFDIAG_H

#include "Kernel.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

/**
 * Jacobian-only kernel on a displacement component: the u <-> eps_bar blocks (SCALAR mean strain,
 * MOOSE's GlobalStrain system), which neither StressDivergenceTensors nor the GlobalStrain scalar
 * kernel assembles. Residual and diagonal Jacobian are zero. Assumes translationally periodic cells.
 */
class GlobalStrainOffDiag : public Kernel
{
public:
  GlobalStrainOffDiag(const InputParameters & parameters);

  static InputParameters validParams();

  virtual void computeResidual() override {}
  virtual void computeJacobian() override {}
  virtual void computeOffDiagJacobian(unsigned int) override {}
  virtual void computeOffDiagJacobianScalar(unsigned int jvar) override;

protected:
  virtual Real computeQpResidual() override { return 0.0; }


private:
  const std::string _base_name;
  const unsigned int _component;
  const unsigned int _global_strain_var;
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;
  DenseMatrix<Number> _ke_copy;
};

#endif //GLOBALSTRAINOFFDIAG_H
