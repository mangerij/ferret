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

#include "GlobalStrainOffDiag.h"
#include "GlobalStrainComponents.h"
#include "MooseVariableScalar.h"
#include "SystemBase.h"
#include "Assembly.h"

registerMooseObject("FerretApp", GlobalStrainOffDiag);

InputParameters
GlobalStrainOffDiag::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Jacobian-only kernel on a displacement component: the u <-> eps_bar "
                             "blocks of MOOSE's GlobalStrain (SCALAR mean strain) system, both "
                             "directions. Use next to StressDivergenceTensors.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "component", "component < 3", "The displacement component this kernel's variable is");
  params.addRequiredCoupledVar("global_strain", "The SCALAR mean strain variable");
  params.addParam<std::string>("base_name", "Material property base name");
  return params;
}

GlobalStrainOffDiag::GlobalStrainOffDiag(const InputParameters & parameters)
  : Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _component(getParam<unsigned int>("component")),
    _global_strain_var(coupledScalar("global_strain")),
    _Jacobian_mult(getMaterialPropertyByName<RankFourTensor>(_base_name + "Jacobian_mult"))
{
}

void
GlobalStrainOffDiag::computeOffDiagJacobianScalar(unsigned int jvar)
{
  // u_c <- eps_bar_b: the mean strain enters sigma like a uniform strain, R_u_c = int sigma_cj dpsi/dx_j,
  //     dR_u_c / d eps_bar_b = int (C:E_b)_cj dpsi/dx_j .
  // eps_bar_b <- u_c: MOOSE's scalar residual is int sigma_(b) for the SINGLE tensor component (b);
  // E_b moves both eps_kl and eps_lk, so (major symmetry of C) the row block is the transpose for
  // the normal components and half of it for the shear components.
  if (jvar != _global_strain_var)
    return;

  MooseVariableScalar & jv = _sys.getScalarVariable(_tid, jvar);

  prepareMatrixTag(_assembly, _var.number(), jvar);
  for (_i = 0; _i < _test.size(); ++_i)
    for (_j = 0; _j < jv.order(); ++_j)
    {
      const auto kl = GlobalStrainComponents::tensorIndices(jv.order(), _j);
      RankTwoTensor E;
      E(kl.first, kl.second) = E(kl.second, kl.first) = 1.0;
      for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      {
        const RankTwoTensor dsigma = _Jacobian_mult[_qp] * E;
        Real r = 0.0;
        for (unsigned int j = 0; j < 3; ++j)
          r += dsigma(_component, j) * _grad_test[_i][_qp](j);
        _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * r;
      }
    }
  accumulateTaggedLocalMatrix();

  _ke_copy = _local_ke;
  prepareMatrixTag(_assembly, jvar, _var.number());
  _ke_copy.get_transpose(_local_ke);
  for (unsigned int b = 0; b < _local_ke.m(); ++b)
    if (GlobalStrainComponents::isShear(jv.order(), b))
      for (unsigned int c = 0; c < _local_ke.n(); ++c)
        _local_ke(b, c) *= 0.5;
  accumulateTaggedLocalMatrix();
}
