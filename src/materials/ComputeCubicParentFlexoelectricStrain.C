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

#include "ComputeCubicParentFlexoelectricStrain.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", ComputeCubicParentFlexoelectricStrain);

/*
  The flexoelectric eigenstrain is eps^F_ij = -A_ijkl dP_k/dx_l with A = s : f, where the
  flexocoupling tensor f is entered through its TENSOR components mu11 = f_1111,
  mu12 = f_1122, mu44 = f_1212 and the compliances in Voigt form (s_1212 = s44/4).
  With this sign, the published f = C : F of Gu et al., PRB 89, 174111 (2014), entered
  directly, gives their tail-to-tail Neel component.  The eigenstrain enters the
  CompositeEigenstrain with weight +1.
  _vals holds TENSOR components in the order xx yy zz yz xz xy (fillFromInputVector).
*/

InputParameters ComputeCubicParentFlexoelectricStrain::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();
  params.addClassDescription("Computes the flexoelectric eigenstrain -(s:f):grad P of a cubic parent phase. It also declares the constant derivative <eigenstrain_name>_dgrad = d eps0/d(grad P), which EigenstrainGradientOffDiag uses for the u<-P Jacobian block.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

ComputeCubicParentFlexoelectricStrain::ComputeCubicParentFlexoelectricStrain(const InputParameters & parameters) :
    ComputeEigenstrainBase(parameters),
   _polar_x_grad(coupledGradient("polar_x")),
   _polar_y_grad(coupledGradient("polar_y")),
   _polar_z_grad(coupledGradient("polar_z")),
   _s11(getMaterialProperty<Real>("s11")),
   _s12(getMaterialProperty<Real>("s12")),
   _s44(getMaterialProperty<Real>("s44")),
   _mu11(getMaterialProperty<Real>("mu11")),
   _mu12(getMaterialProperty<Real>("mu12")),
   _mu44(getMaterialProperty<Real>("mu44")),
   _vals(6),
   _flexoelectric_strain(),
   _deigenstrain_dgrad(declareProperty<RankFourTensor>(_eigenstrain_name + "_dgrad"))
{
}

void
ComputeCubicParentFlexoelectricStrain::computeQpEigenstrain()
{
  RankTwoTensor grad_p;
  for (unsigned int l = 0; l < 3; ++l)
  {
    grad_p(0, l) = _polar_x_grad[_qp](l);
    grad_p(1, l) = _polar_y_grad[_qp](l);
    grad_p(2, l) = _polar_z_grad[_qp](l);
  }
  _eigenstrain[_qp] = flexoStrain(grad_p);

  // eps0 is linear in grad P, so its derivative is eps0 evaluated on the unit gradients
  for (unsigned int k = 0; k < 3; ++k)
    for (unsigned int l = 0; l < 3; ++l)
    {
      RankTwoTensor unit;
      unit(k, l) = 1.0;
      const RankTwoTensor d = flexoStrain(unit);
      for (unsigned int m = 0; m < 3; ++m)
        for (unsigned int n = 0; n < 3; ++n)
          _deigenstrain_dgrad[_qp](m, n, k, l) = d(m, n);
    }
}

RankTwoTensor
ComputeCubicParentFlexoelectricStrain::flexoStrain(const RankTwoTensor & grad_p)
{
  _vals[0] = -grad_p(0, 0)*(_mu11[_qp]*_s11[_qp] + 2*_mu12[_qp]*_s12[_qp]) - (grad_p(1, 1) + grad_p(2, 2))*(_mu11[_qp]*_s12[_qp] + _mu12[_qp]*_s11[_qp]
      + _mu12[_qp]*_s12[_qp]);  // eps_1   11
  _vals[1] = -grad_p(1, 1)*(_mu11[_qp]*_s11[_qp] + 2*_mu12[_qp]*_s12[_qp]) - (grad_p(0, 0) + grad_p(2, 2))*(_mu11[_qp]*_s12[_qp] + _mu12[_qp]*_s11[_qp]
      + _mu12[_qp]*_s12[_qp]);  // eps_2   22
  _vals[2] = -grad_p(2, 2)*(_mu11[_qp]*_s11[_qp] + 2*_mu12[_qp]*_s12[_qp]) - (grad_p(0, 0) + grad_p(1, 1))*(_mu11[_qp]*_s12[_qp] + _mu12[_qp]*_s11[_qp]
      + _mu12[_qp]*_s12[_qp]);  // eps_3   33
  _vals[3] = -0.5*_mu44[_qp]*_s44[_qp]*(grad_p(1, 2) + grad_p(2, 1));  // eps_4   23
  _vals[4] = -0.5*_mu44[_qp]*_s44[_qp]*(grad_p(0, 2) + grad_p(2, 0));  // eps_5   13
  _vals[5] = -0.5*_mu44[_qp]*_s44[_qp]*(grad_p(0, 1) + grad_p(1, 0));  // eps_6   12
  _flexoelectric_strain.fillFromInputVector(_vals);
  return _flexoelectric_strain;
}
