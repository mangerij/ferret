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

#ifndef EIGENSTRAINGRADIENTOFFDIAG_H
#define EIGENSTRAINGRADIENTOFFDIAG_H

#include "Kernel.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

/**
 * Jacobian-only kernel on a displacement component: the u <- eta block of an eigenstrain that is
 * linear in the GRADIENT of an order parameter eta (e.g. the flexoelectric strain, eps0 = G grad P).
 * StressDivergenceTensors only differentiates eigenstrains w.r.t. variable values, so it cannot
 * supply this block. Residual and diagonal Jacobian are zero; StressDivergenceTensors keeps both.
 */
class EigenstrainGradientOffDiag : public Kernel
{
public:
  EigenstrainGradientOffDiag(const InputParameters & parameters);

  static InputParameters validParams();

  virtual void computeResidual() override {}
  virtual void computeJacobian() override {}

protected:
  virtual Real computeQpResidual() override { return 0.0; }
  virtual void precalculateOffDiagJacobian(unsigned int jvar) override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  /// which order-parameter component jvar is, or libMesh::invalid_uint
  unsigned int orderParameterComponent(unsigned int jvar) const;

  const std::string _base_name;
  const unsigned int _component;
  const unsigned int _n_eta;
  std::vector<unsigned int> _eta_var;
  const Real _prefactor;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  /// d eps0_mn / d(d eta_k / dx_l), declared by the eigenstrain material as <eigenstrain_name>_dgrad
  const MaterialProperty<RankFourTensor> & _deigenstrain_dgrad;
  /// A_jl = C_(component) j m n G_mnkl for the column component k, per qp
  std::vector<RankTwoTensor> _A;
};

#endif //EIGENSTRAINGRADIENTOFFDIAG_H
