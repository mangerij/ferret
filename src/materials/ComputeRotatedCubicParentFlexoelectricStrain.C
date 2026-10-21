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

#include "ComputeRotatedCubicParentFlexoelectricStrain.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", ComputeRotatedCubicParentFlexoelectricStrain);

/*
  The flexoelectric eigenstrain is eps^F_ij = -A_ijkl dP_k/dx_l with A = s : f, where the
  flexocoupling tensor f is entered through its TENSOR components mu11 = f_1111,
  mu12 = f_1122, mu44 = f_1212 and the compliances in Voigt form (s_1212 = s44/4).
  With this sign, the published f = C : F of Gu et al., PRB 89, 174111 (2014), entered
  directly, gives their tail-to-tail Neel component.  The eigenstrain enters the
  CompositeEigenstrain with weight +1.
  _vals holds TENSOR components in the order xx yy zz yz xz xy (fillFromInputVector).
*/

InputParameters ComputeRotatedCubicParentFlexoelectricStrain::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();
  params.addClassDescription("Computes the flexoelectric eigenstrain -(s:f):grad P of a cubic parent phase whose crystal axes are turned by the material property twist about z.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

ComputeRotatedCubicParentFlexoelectricStrain::ComputeRotatedCubicParentFlexoelectricStrain(const InputParameters & parameters) :
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
   _twist(getMaterialProperty<Real>("twist")),
   _vals(6),
   _flexoelectric_strain()
{
}

void
ComputeRotatedCubicParentFlexoelectricStrain::computeQpEigenstrain()
{
  const Real ct = std::cos(_twist[_qp]);
  const Real st = std::sin(_twist[_qp]);
  const Real g00 = _polar_x_grad[_qp](0)*Utility::pow<2>(ct) + _polar_x_grad[_qp](1)*ct*st + _polar_y_grad[_qp](0)*ct*st + _polar_y_grad[_qp](1)*Utility::pow<2>(st);
  const Real g01 = -_polar_x_grad[_qp](0)*ct*st + _polar_x_grad[_qp](1)*Utility::pow<2>(ct) - _polar_y_grad[_qp](0)*Utility::pow<2>(st) + _polar_y_grad[_qp](1)*ct*st;
  const Real g02 = _polar_x_grad[_qp](2)*ct + _polar_y_grad[_qp](2)*st;
  const Real g10 = -_polar_x_grad[_qp](0)*ct*st - _polar_x_grad[_qp](1)*Utility::pow<2>(st) + _polar_y_grad[_qp](0)*Utility::pow<2>(ct) + _polar_y_grad[_qp](1)*ct*st;
  const Real g11 = _polar_x_grad[_qp](0)*Utility::pow<2>(st) - _polar_x_grad[_qp](1)*ct*st - _polar_y_grad[_qp](0)*ct*st + _polar_y_grad[_qp](1)*Utility::pow<2>(ct);
  const Real g12 = -_polar_x_grad[_qp](2)*st + _polar_y_grad[_qp](2)*ct;
  const Real g20 = _polar_z_grad[_qp](0)*ct + _polar_z_grad[_qp](1)*st;
  const Real g21 = -_polar_z_grad[_qp](0)*st + _polar_z_grad[_qp](1)*ct;
  const Real g22 = _polar_z_grad[_qp](2);
  // the eigenstrain in the crystal frame, then turned back: eps^F = R eps^F_c R^T
  const Real f00 = -g00*(_mu11[_qp]*_s11[_qp] + 2*_mu12[_qp]*_s12[_qp]) - (g11 + g22)*(_mu11[_qp]*_s12[_qp] + _mu12[_qp]*_s11[_qp] + _mu12[_qp]*_s12[_qp]);
  const Real f01 = -0.5*_mu44[_qp]*_s44[_qp]*(g01 + g10);
  const Real f02 = -0.5*_mu44[_qp]*_s44[_qp]*(g02 + g20);
  const Real f11 = -g11*(_mu11[_qp]*_s11[_qp] + 2*_mu12[_qp]*_s12[_qp]) - (g00 + g22)*(_mu11[_qp]*_s12[_qp] + _mu12[_qp]*_s11[_qp] + _mu12[_qp]*_s12[_qp]);
  const Real f12 = -0.5*_mu44[_qp]*_s44[_qp]*(g12 + g21);
  const Real f22 = -g22*(_mu11[_qp]*_s11[_qp] + 2*_mu12[_qp]*_s12[_qp]) - (g00 + g11)*(_mu11[_qp]*_s12[_qp] + _mu12[_qp]*_s11[_qp] + _mu12[_qp]*_s12[_qp]);

  _vals[0] = Utility::pow<2>(ct)*f00 - 2*ct*f01*st + f11*Utility::pow<2>(st);  // eps_1   11
  _vals[1] = Utility::pow<2>(ct)*f11 + 2*ct*f01*st + f00*Utility::pow<2>(st);  // eps_2   22
  _vals[2] = f22;  // eps_3   33
  _vals[3] = ct*f12 + f02*st;  // eps_4   23
  _vals[4] = ct*f02 - f12*st;  // eps_5   13
  _vals[5] = Utility::pow<2>(ct)*f01 + ct*f00*st - ct*f11*st - f01*Utility::pow<2>(st);  // eps_6   12
  _flexoelectric_strain.fillFromInputVector(_vals);
  _eigenstrain[_qp] = _flexoelectric_strain;
}
