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

#include "EigenstrainGradientOffDiag.h"

registerMooseObject("FerretApp", EigenstrainGradientOffDiag);

InputParameters
EigenstrainGradientOffDiag::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "Jacobian-only kernel on a displacement component: the u<-eta block of an eigenstrain that is "
      "linear in grad eta (e.g. the flexoelectric strain), read from the material property "
      "<eigenstrain_name>_dgrad. Use next to StressDivergenceTensors, which supplies the residual.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "component", "component < 3", "The displacement component this kernel's variable is");
  params.addRequiredCoupledVar(
      "order_parameter",
      "The order-parameter components eta_k, in the order of the k index of <eigenstrain_name>_dgrad");
  params.addRequiredParam<std::string>(
      "eigenstrain_name", "The eigenstrain whose gradient derivative <eigenstrain_name>_dgrad is read");
  params.addParam<Real>("prefactor",
                        1.0,
                        "Constant weight of this eigenstrain in the total eigenstrain (e.g. its "
                        "CompositeEigenstrain weight)");
  params.addParam<std::string>("base_name", "Material property base name");
  return params;
}

EigenstrainGradientOffDiag::EigenstrainGradientOffDiag(const InputParameters & parameters)
  : Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _component(getParam<unsigned int>("component")),
    _n_eta(coupledComponents("order_parameter")),
    _eta_var(_n_eta),
    _prefactor(getParam<Real>("prefactor")),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")),
    _deigenstrain_dgrad(getMaterialPropertyByName<RankFourTensor>(
        _base_name + getParam<std::string>("eigenstrain_name") + "_dgrad"))
{
  if (_n_eta > 3)
    paramError("order_parameter", "at most three components");
  for (unsigned int k = 0; k < _n_eta; ++k)
    _eta_var[k] = coupled("order_parameter", k);
}

unsigned int
EigenstrainGradientOffDiag::orderParameterComponent(unsigned int jvar) const
{
  for (unsigned int k = 0; k < _n_eta; ++k)
    if (jvar == _eta_var[k])
      return k;
  return libMesh::invalid_uint;
}

void
EigenstrainGradientOffDiag::precalculateOffDiagJacobian(unsigned int jvar)
{
  // R_u_i = int sigma_ij dpsi/dx_j,  sigma = C : (eps - prefactor G grad eta), so
  //   dR_u_i / d eta_k [phi] = -prefactor int dpsi/dx_j C_ijmn G_mnkl dphi/dx_l .
  // A = C_i.mn G_mnkl depends on the qp only: contract once here, not inside the i,j loops.
  const unsigned int k = orderParameterComponent(jvar);
  if (k == libMesh::invalid_uint)
    return;

  _A.resize(_qrule->n_points());
  for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
  {
    _A[qp].zero();
    for (unsigned int j = 0; j < 3; ++j)
      for (unsigned int l = 0; l < 3; ++l)
        for (unsigned int m = 0; m < 3; ++m)
          for (unsigned int n = 0; n < 3; ++n)
            _A[qp](j, l) +=
                _elasticity_tensor[qp](_component, j, m, n) * _deigenstrain_dgrad[qp](m, n, k, l);
  }
}

Real
EigenstrainGradientOffDiag::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (orderParameterComponent(jvar) == libMesh::invalid_uint)
    return 0.0;

  return -_prefactor * (_grad_test[_i][_qp] * (_A[_qp] * _grad_phi[_j][_qp]));
}
