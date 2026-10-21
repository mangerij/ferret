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

#include "CubicParentFlexoelectricFullEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", CubicParentFlexoelectricFullEnergy);

/*
  The energy whose P-variation is CubicParentFlexoelectricFullPDerivative.
  The flexoelectric eigenstrain is eps^F_ij = -A_ijkl dP_k/dx_l with A = s : f, where the
  flexocoupling tensor f is entered through its TENSOR components mu11 = f_1111,
  mu12 = f_1122, mu44 = f_1212 and the compliances in Voigt form (s_1212 = s44/4).
  With this sign, the published f = C : F of Gu et al., PRB 89, 174111 (2014), entered
  directly, gives their tail-to-tail Neel component.  The eigenstrain enters the
  CompositeEigenstrain with weight +1.
*/

InputParameters CubicParentFlexoelectricFullEnergy::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates the volume integral of the full elastic energy density 1/2 (eps - eps^Q - eps^F):C:(eps - eps^Q - eps^F) of a cubic parent phase, electrostrictive and flexoelectric eigenstrains included.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

CubicParentFlexoelectricFullEnergy::CubicParentFlexoelectricFullEnergy(const InputParameters & parameters) :
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
   _energy_scale(getParam<Real>("energy_scale")),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
   _strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_strain"))
{
}

Real
CubicParentFlexoelectricFullEnergy::computeQpIntegral()
{
  return _energy_scale*((0.5)*_C11[_qp]*(Utility::pow<2>(_strain[_qp](0,0)) + Utility::pow<2>(_strain[_qp](1,1)) + Utility::pow<2>(_strain[_qp](2,2))) + _C12[_qp]*(_strain[_qp](0,0)*_strain[_qp](1,1)
      + _strain[_qp](0,0)*_strain[_qp](2,2) + _strain[_qp](1,1)*_strain[_qp](2,2)) - 2*_C44[_qp]*_Q44[_qp]*_mu44[_qp]*_s44[_qp]*(_polar_x_grad[_qp](1)*_polar_x[_qp]*_polar_y[_qp]
      + _polar_x_grad[_qp](2)*_polar_x[_qp]*_polar_z[_qp] + _polar_y_grad[_qp](0)*_polar_x[_qp]*_polar_y[_qp] + _polar_y_grad[_qp](2)*_polar_y[_qp]*_polar_z[_qp]
      + _polar_z_grad[_qp](0)*_polar_x[_qp]*_polar_z[_qp] + _polar_z_grad[_qp](1)*_polar_y[_qp]*_polar_z[_qp]) - 4*_C44[_qp]*_Q44[_qp]*(_strain[_qp](0,1)*_polar_x[_qp]*_polar_y[_qp]
      + _strain[_qp](0,2)*_polar_x[_qp]*_polar_z[_qp] + _strain[_qp](1,2)*_polar_y[_qp]*_polar_z[_qp]) + _C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(_polar_x_grad[_qp](1)*_polar_y_grad[_qp](0)
      + _polar_x_grad[_qp](2)*_polar_z_grad[_qp](0) + _polar_y_grad[_qp](2)*_polar_z_grad[_qp](1)) + (0.5)*_C44[_qp]*Utility::pow<2>(_mu44[_qp])*Utility::pow<2>(_s44[_qp])*(Utility::pow<2>(_polar_x_grad[_qp](1))
      + Utility::pow<2>(_polar_x_grad[_qp](2)) + Utility::pow<2>(_polar_y_grad[_qp](0)) + Utility::pow<2>(_polar_y_grad[_qp](2)) + Utility::pow<2>(_polar_z_grad[_qp](0))
      + Utility::pow<2>(_polar_z_grad[_qp](1))) + 2*_C44[_qp]*_mu44[_qp]*_s44[_qp]*(_strain[_qp](0,1)*_polar_x_grad[_qp](1) + _strain[_qp](0,1)*_polar_y_grad[_qp](0)
      + _strain[_qp](0,2)*_polar_x_grad[_qp](2) + _strain[_qp](0,2)*_polar_z_grad[_qp](0) + _strain[_qp](1,2)*_polar_y_grad[_qp](2) + _strain[_qp](1,2)*_polar_z_grad[_qp](1))
      + 2*_C44[_qp]*(Utility::pow<2>(_strain[_qp](0,1)) + Utility::pow<2>(_strain[_qp](0,2)) + Utility::pow<2>(_strain[_qp](1,2))) - (_C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp])*(_strain[_qp](0,0)*Utility::pow<2>(_polar_x[_qp])
      + _strain[_qp](1,1)*Utility::pow<2>(_polar_y[_qp]) + _strain[_qp](2,2)*Utility::pow<2>(_polar_z[_qp])) + (Utility::pow<2>(_polar_x_grad[_qp](0)) + Utility::pow<2>(_polar_y_grad[_qp](1))
      + Utility::pow<2>(_polar_z_grad[_qp](2)))*((0.5)*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp])
      + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp] + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp])
      + 2*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp] + 3*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 6*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 5*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) + (Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_z[_qp]))*((0.5)*_C11[_qp]*Utility::pow<2>(_Q11[_qp])
      + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + _C12[_qp]*Utility::pow<2>(_Q12[_qp])) - (_C11[_qp]*_Q12[_qp] + _C12[_qp]*_Q11[_qp]
      + _C12[_qp]*_Q12[_qp])*(_strain[_qp](0,0)*Utility::pow<2>(_polar_y[_qp]) + _strain[_qp](0,0)*Utility::pow<2>(_polar_z[_qp]) + _strain[_qp](1,1)*Utility::pow<2>(_polar_x[_qp])
      + _strain[_qp](1,1)*Utility::pow<2>(_polar_z[_qp]) + _strain[_qp](2,2)*Utility::pow<2>(_polar_x[_qp]) + _strain[_qp](2,2)*Utility::pow<2>(_polar_y[_qp]))
      + (_strain[_qp](0,0)*_polar_x_grad[_qp](0) + _strain[_qp](1,1)*_polar_y_grad[_qp](1) + _strain[_qp](2,2)*_polar_z_grad[_qp](2))*(_C11[_qp]*_mu11[_qp]*_s11[_qp]
      + 2*_C11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C12[_qp]*_mu12[_qp]*_s12[_qp]) + (_polar_x_grad[_qp](0)*_polar_y_grad[_qp](1)
      + _polar_x_grad[_qp](0)*_polar_z_grad[_qp](2) + _polar_y_grad[_qp](1)*_polar_z_grad[_qp](2))*(2*_C11[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + _C11[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 4*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 6*_C11[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + _C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 6*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 5*_C11[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp]) + _C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s11[_qp]) + 2*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*_s11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*Utility::pow<2>(_mu11[_qp])*Utility::pow<2>(_s12[_qp]) + 2*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s11[_qp]) + 12*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*_s11[_qp]*_s12[_qp]
      + 10*_C12[_qp]*_mu11[_qp]*_mu12[_qp]*Utility::pow<2>(_s12[_qp]) + 3*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s11[_qp]) + 10*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*_s11[_qp]*_s12[_qp]
      + 11*_C12[_qp]*Utility::pow<2>(_mu12[_qp])*Utility::pow<2>(_s12[_qp])) - (_polar_x_grad[_qp](0)*Utility::pow<2>(_polar_x[_qp]) + _polar_y_grad[_qp](1)*Utility::pow<2>(_polar_y[_qp])
      + _polar_z_grad[_qp](2)*Utility::pow<2>(_polar_z[_qp]))*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 2*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp]
      + 2*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp] + 2*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp]
      + 6*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]) + (Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp])
      + Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*(2*_C11[_qp]*_Q11[_qp]*_Q12[_qp] + _C11[_qp]*Utility::pow<2>(_Q12[_qp]) + _C12[_qp]*Utility::pow<2>(_Q11[_qp])
      + 2*_C12[_qp]*_Q11[_qp]*_Q12[_qp] + 3*_C12[_qp]*Utility::pow<2>(_Q12[_qp]) + 2*_C44[_qp]*Utility::pow<2>(_Q44[_qp])) + (_strain[_qp](0,0)*_polar_y_grad[_qp](1)
      + _strain[_qp](0,0)*_polar_z_grad[_qp](2) + _strain[_qp](1,1)*_polar_x_grad[_qp](0) + _strain[_qp](1,1)*_polar_z_grad[_qp](2) + _strain[_qp](2,2)*_polar_x_grad[_qp](0)
      + _strain[_qp](2,2)*_polar_y_grad[_qp](1))*(_C11[_qp]*_mu11[_qp]*_s12[_qp] + _C11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_mu11[_qp]*_s11[_qp]
      + _C12[_qp]*_mu11[_qp]*_s12[_qp] + _C12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_mu12[_qp]*_s12[_qp]) - (_polar_x_grad[_qp](0)*Utility::pow<2>(_polar_y[_qp])
      + _polar_x_grad[_qp](0)*Utility::pow<2>(_polar_z[_qp]) + _polar_y_grad[_qp](1)*Utility::pow<2>(_polar_x[_qp]) + _polar_y_grad[_qp](1)*Utility::pow<2>(_polar_z[_qp])
      + _polar_z_grad[_qp](2)*Utility::pow<2>(_polar_x[_qp]) + _polar_z_grad[_qp](2)*Utility::pow<2>(_polar_y[_qp]))*(_C11[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + _C11[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + _C11[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + _C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C11[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s11[_qp] + _C12[_qp]*_Q11[_qp]*_mu11[_qp]*_s12[_qp]
      + _C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q11[_qp]*_mu12[_qp]*_s12[_qp] + _C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s11[_qp] + 3*_C12[_qp]*_Q12[_qp]*_mu11[_qp]*_s12[_qp]
      + 3*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s11[_qp] + 5*_C12[_qp]*_Q12[_qp]*_mu12[_qp]*_s12[_qp]));
}
