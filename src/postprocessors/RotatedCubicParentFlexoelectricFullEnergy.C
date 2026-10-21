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

#include "RotatedCubicParentFlexoelectricFullEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", RotatedCubicParentFlexoelectricFullEnergy);

/*
  The energy whose P-variation is RotatedCubicParentFlexoelectricFullPDerivative.
  The flexoelectric eigenstrain is eps^F_ij = -A_ijkl dP_k/dx_l with A = s : f, where the
  flexocoupling tensor f is entered through its TENSOR components mu11 = f_1111,
  mu12 = f_1122, mu44 = f_1212 and the compliances in Voigt form (s_1212 = s44/4).
  With this sign, the published f = C : F of Gu et al., PRB 89, 174111 (2014), entered
  directly, gives their tail-to-tail Neel component.  The eigenstrain enters the
  CompositeEigenstrain with weight +1.
*/

InputParameters RotatedCubicParentFlexoelectricFullEnergy::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates the volume integral of the full elastic energy density 1/2 (eps - eps^Q - eps^F):C:(eps - eps^Q - eps^F) of a cubic parent phase, electrostrictive and flexoelectric eigenstrains included. The crystal axes are turned by the material property twist about z.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

RotatedCubicParentFlexoelectricFullEnergy::RotatedCubicParentFlexoelectricFullEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _polar_x_grad(coupledGradient("polar_x")),
   _polar_y_grad(coupledGradient("polar_y")),
   _polar_z_grad(coupledGradient("polar_z")),
   _C11(getMaterialProperty<Real>("C11")),
   _C12(getMaterialProperty<Real>("C12")),
   _C44(getMaterialProperty<Real>("C44")),
   _Q11(getMaterialProperty<Real>("Q11")),
   _Q12(getMaterialProperty<Real>("Q12")),
   _Q44(getMaterialProperty<Real>("Q44")),
   _s11(getMaterialProperty<Real>("s11")),
   _s12(getMaterialProperty<Real>("s12")),
   _s44(getMaterialProperty<Real>("s44")),
   _mu11(getMaterialProperty<Real>("mu11")),
   _mu12(getMaterialProperty<Real>("mu12")),
   _mu44(getMaterialProperty<Real>("mu44")),
   _twist(getMaterialProperty<Real>("twist")),
   _energy_scale(getParam<Real>("energy_scale")),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
   _strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_strain"))
{
}

Real
RotatedCubicParentFlexoelectricFullEnergy::computeQpIntegral()
{
  const Real ct = std::cos(_twist[_qp]);
  const Real st = std::sin(_twist[_qp]);
  // polarization, its gradient and the total strain in the crystal frame
  const Real p0 = _polar_x[_qp]*ct + _polar_y[_qp]*st;
  const Real p1 = -_polar_x[_qp]*st + _polar_y[_qp]*ct;
  const Real p2 = _polar_z[_qp];
  const Real g00 = _polar_x_grad[_qp](0)*Utility::pow<2>(ct) + _polar_x_grad[_qp](1)*ct*st + _polar_y_grad[_qp](0)*ct*st + _polar_y_grad[_qp](1)*Utility::pow<2>(st);
  const Real g01 = -_polar_x_grad[_qp](0)*ct*st + _polar_x_grad[_qp](1)*Utility::pow<2>(ct) - _polar_y_grad[_qp](0)*Utility::pow<2>(st) + _polar_y_grad[_qp](1)*ct*st;
  const Real g02 = _polar_x_grad[_qp](2)*ct + _polar_y_grad[_qp](2)*st;
  const Real g10 = -_polar_x_grad[_qp](0)*ct*st - _polar_x_grad[_qp](1)*Utility::pow<2>(st) + _polar_y_grad[_qp](0)*Utility::pow<2>(ct) + _polar_y_grad[_qp](1)*ct*st;
  const Real g11 = _polar_x_grad[_qp](0)*Utility::pow<2>(st) - _polar_x_grad[_qp](1)*ct*st - _polar_y_grad[_qp](0)*ct*st + _polar_y_grad[_qp](1)*Utility::pow<2>(ct);
  const Real g12 = -_polar_x_grad[_qp](2)*st + _polar_y_grad[_qp](2)*ct;
  const Real g20 = _polar_z_grad[_qp](0)*ct + _polar_z_grad[_qp](1)*st;
  const Real g21 = -_polar_z_grad[_qp](0)*st + _polar_z_grad[_qp](1)*ct;
  const Real g22 = _polar_z_grad[_qp](2);
  const Real e00 = _strain[_qp](0,0)*Utility::pow<2>(ct) + 2*_strain[_qp](0,1)*ct*st + _strain[_qp](1,1)*Utility::pow<2>(st);
  const Real e01 = -_strain[_qp](0,0)*ct*st + _strain[_qp](0,1)*Utility::pow<2>(ct) - _strain[_qp](0,1)*Utility::pow<2>(st) + _strain[_qp](1,1)*ct*st;
  const Real e02 = _strain[_qp](0,2)*ct + _strain[_qp](1,2)*st;
  const Real e11 = _strain[_qp](0,0)*Utility::pow<2>(st) - 2*_strain[_qp](0,1)*ct*st + _strain[_qp](1,1)*Utility::pow<2>(ct);
  const Real e12 = -_strain[_qp](0,2)*st + _strain[_qp](1,2)*ct;
  const Real e22 = _strain[_qp](2,2);

  return _energy_scale*((0.5)*_C11[_qp]*(Utility::pow<2>(e00) + Utility::pow<2>(e11) + Utility::pow<2>(e22)) + _C12[_qp]*(e00*e11 + e00*e22 + e11*e22) - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(g01*p0*p1
      + g02*p0*p2 + g10*p0*p1 + g12*p1*p2 + g20*p0*p2 + g21*p1*p2) - 4*_C44[_qp]*_Q44[_qp]*(e01*p0*p1 + e02*p0*p2 + e12*p1*p2) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(g01*g10
      + g02*g20 + g12*g21) + (0.5)*_C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(Utility::pow<2>(g01) + Utility::pow<2>(g02) + Utility::pow<2>(g10)
      + Utility::pow<2>(g12) + Utility::pow<2>(g20) + Utility::pow<2>(g21)) + 2*_C44[_qp]*_mu44[_qp]*_s44[_qp]*(e01*g01 + e01*g10 + e02*g02 + e02*g20 + e12*g12
      + e12*g21) + 2*_C44[_qp]*(Utility::pow<2>(e01) + Utility::pow<2>(e02) + Utility::pow<2>(e12)) - (_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp])*(e00*Utility::pow<2>(p0)
      + e11*Utility::pow<2>(p1) + e22*Utility::pow<2>(p2)) + (Utility::pow<2>(g00) + Utility::pow<2>(g11) + Utility::pow<2>(g22))*((0.5)*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp])
      + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 3*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])
      + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 6*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) + (Utility::pow<4>(p0)
      + Utility::pow<4>(p1) + Utility::pow<4>(p2))*((0.5)*_C11[_qp]*Utility::pow<2>(_Q11[_qp]) + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C12[_qp]*Utility::pow<2>(_Q12[_qp])) - (_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp] + _C12[_qp]*_Q12[_qp])*(e00*Utility::pow<2>(p1) + e00*Utility::pow<2>(p2)
      + e11*Utility::pow<2>(p0) + e11*Utility::pow<2>(p2) + e22*Utility::pow<2>(p0) + e22*Utility::pow<2>(p1)) + (e00*g00 + e11*g11 + e22*g22)*(_C11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]) + (g00*g11 + g00*g22
      + g11*g22)*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp])
      + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - (g00*Utility::pow<2>(p0)
      + g11*Utility::pow<2>(p1) + g22*Utility::pow<2>(p2))*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp]
      + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + (Utility::pow<2>(p0)*Utility::pow<2>(p1) + Utility::pow<2>(p0)*Utility::pow<2>(p2) + Utility::pow<2>(p1)*Utility::pow<2>(p2))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp]
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp])
      + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + (e00*g11 + e00*g22 + e11*g00 + e11*g22 + e22*g00 + e22*g11)*(_C11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_mu12[_qp]*_s11[_qp]
      + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp])
      - (g00*Utility::pow<2>(p1) + g00*Utility::pow<2>(p2) + g11*Utility::pow<2>(p0) + g11*Utility::pow<2>(p2) + g22*Utility::pow<2>(p0) + g22*Utility::pow<2>(p1))*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]));
}
