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

#include "ComputeRotatedCubicParentElectrostrictiveStrain.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ComputeRotatedCubicParentElectrostrictiveStrain);

InputParameters ComputeRotatedCubicParentElectrostrictiveStrain::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();
  params.addClassDescription("Computes the electrostrictive eigenstrain eps_ij = Q_ijkl P_k P_l of a cubic parent phase whose crystal axes are turned by the material property twist about z. It also declares its derivatives w.r.t. the coupled polarization components, so StressDivergenceTensors assembles the u<-P Jacobian block when given eigenstrain_names = <this eigenstrain_name> and coupled_variables = <the polarization components>.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

ComputeRotatedCubicParentElectrostrictiveStrain::ComputeRotatedCubicParentElectrostrictiveStrain(const InputParameters & parameters) :
    DerivativeMaterialInterface<ComputeEigenstrainBase>(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _Q11(getMaterialProperty<Real>("Q11")),
  _Q12(getMaterialProperty<Real>("Q12")),
  _Q44(getMaterialProperty<Real>("Q44")),
  _twist(getMaterialProperty<Real>("twist")),
  _vals(6),
  _rotated_cubic_parent_electrostrictive_strain(),
  _deigenstrain(3, nullptr)
{
  const std::vector<std::string> comps = {"polar_x", "polar_y", "polar_z"};
  for (unsigned int k = 0; k < 3; ++k)
    if (isCoupled(comps[k]))
      _deigenstrain[k] = &declarePropertyDerivative<RankTwoTensor>(_eigenstrain_name, coupledName(comps[k]));
}

void
ComputeRotatedCubicParentElectrostrictiveStrain::computeQpEigenstrain()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  _vals[0] = ((3*Utility::pow<2>(_polar_x[_qp])*_Q11[_qp] + Utility::pow<2>(_polar_y[_qp])*_Q11[_qp] + Utility::pow<2>(_polar_x[_qp])*_Q12[_qp] + 3*Utility::pow<2>(_polar_y[_qp])*_Q12[_qp] + 4*Utility::pow<2>(_polar_z[_qp])*_Q12[_qp] + Utility::pow<2>(_polar_x[_qp])*_Q44[_qp] - Utility::pow<2>(_polar_y[_qp])*_Q44[_qp] + 
     (Utility::pow<2>(_polar_x[_qp]) - Utility::pow<2>(_polar_y[_qp]))*(_Q11[_qp] - _Q12[_qp] - _Q44[_qp])*std::cos(4*_twist[_qp]) + 2*_polar_x[_qp]*_polar_y[_qp]*(_Q11[_qp] - _Q12[_qp] - _Q44[_qp])*std::sin(4*_twist[_qp]))/4.);  // eps_1
  _vals[1] = ((Utility::pow<2>(_polar_x[_qp])*_Q11[_qp] + 3*Utility::pow<2>(_polar_y[_qp])*_Q11[_qp] + 3*Utility::pow<2>(_polar_x[_qp])*_Q12[_qp] + Utility::pow<2>(_polar_y[_qp])*_Q12[_qp] + 4*Utility::pow<2>(_polar_z[_qp])*_Q12[_qp] - Utility::pow<2>(_polar_x[_qp])*_Q44[_qp] + Utility::pow<2>(_polar_y[_qp])*_Q44[_qp] - 
     (Utility::pow<2>(_polar_x[_qp]) - Utility::pow<2>(_polar_y[_qp]))*(_Q11[_qp] - _Q12[_qp] - _Q44[_qp])*std::cos(4*_twist[_qp]) - 2*_polar_x[_qp]*_polar_y[_qp]*(_Q11[_qp] - _Q12[_qp] - _Q44[_qp])*std::sin(4*_twist[_qp]))/4.);  // eps_2
  _vals[2] = (Utility::pow<2>(_polar_z[_qp])*_Q11[_qp] + (Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*_Q12[_qp]);  // eps_3
  _vals[3] = (_polar_y[_qp]*_polar_z[_qp]*_Q44[_qp]);  // eps_4   23
  _vals[4] = (_polar_x[_qp]*_polar_z[_qp]*_Q44[_qp]);  // eps_5   13
  _vals[5] = ((2*_polar_x[_qp]*_polar_y[_qp]*(_Q11[_qp] - _Q12[_qp] + _Q44[_qp]) - 2*_polar_x[_qp]*_polar_y[_qp]*(_Q11[_qp] - _Q12[_qp] - _Q44[_qp])*std::cos(4*_twist[_qp]) + (Utility::pow<2>(_polar_x[_qp]) - Utility::pow<2>(_polar_y[_qp]))*(_Q11[_qp] - _Q12[_qp] - _Q44[_qp])*std::sin(4*_twist[_qp]))/4.);  // eps_6   12
  _rotated_cubic_parent_electrostrictive_strain.fillFromInputVector(_vals);
  _eigenstrain[_qp] = _rotated_cubic_parent_electrostrictive_strain;

  // d eps0 / d P_k, differentiated (sympy) from the expressions above
  if (_deigenstrain[0])
  {
    _vals[0] = (1.0/2.0)*(_Q11[_qp]*_polar_x[_qp]*std::cos(4*_twist[_qp]) + 3*_Q11[_qp]*_polar_x[_qp] + _Q11[_qp]*_polar_y[_qp]*std::sin(4*_twist[_qp]) - _Q12[_qp]*_polar_x[_qp]*std::cos(4*_twist[_qp]) + _Q12[_qp]*_polar_x[_qp] - _Q12[_qp]*_polar_y[_qp]*std::sin(4*_twist[_qp]) - _Q44[_qp]*_polar_x[_qp]*std::cos(4*_twist[_qp]) + _Q44[_qp]*_polar_x[_qp] - _Q44[_qp]*_polar_y[_qp]*std::sin(4*_twist[_qp]));  // d eps_1
    _vals[1] = (1.0/2.0)*(-_Q11[_qp]*_polar_x[_qp]*std::cos(4*_twist[_qp]) + _Q11[_qp]*_polar_x[_qp] - _Q11[_qp]*_polar_y[_qp]*std::sin(4*_twist[_qp]) + _Q12[_qp]*_polar_x[_qp]*std::cos(4*_twist[_qp]) + 3*_Q12[_qp]*_polar_x[_qp] + _Q12[_qp]*_polar_y[_qp]*std::sin(4*_twist[_qp]) + _Q44[_qp]*_polar_x[_qp]*std::cos(4*_twist[_qp]) - _Q44[_qp]*_polar_x[_qp] + _Q44[_qp]*_polar_y[_qp]*std::sin(4*_twist[_qp]));  // d eps_2
    _vals[2] = 2*_Q12[_qp]*_polar_x[_qp];  // d eps_3
    _vals[3] = 0;  // d eps_4   23
    _vals[4] = _Q44[_qp]*_polar_z[_qp];  // d eps_5   13
    _vals[5] = (1.0/2.0)*(_Q11[_qp]*_polar_x[_qp]*std::sin(4*_twist[_qp]) - _Q11[_qp]*_polar_y[_qp]*std::cos(4*_twist[_qp]) + _Q11[_qp]*_polar_y[_qp] - _Q12[_qp]*_polar_x[_qp]*std::sin(4*_twist[_qp]) + _Q12[_qp]*_polar_y[_qp]*std::cos(4*_twist[_qp]) - _Q12[_qp]*_polar_y[_qp] - _Q44[_qp]*_polar_x[_qp]*std::sin(4*_twist[_qp]) + _Q44[_qp]*_polar_y[_qp]*std::cos(4*_twist[_qp]) + _Q44[_qp]*_polar_y[_qp]);  // d eps_6   12
    _rotated_cubic_parent_electrostrictive_strain.fillFromInputVector(_vals);
    (*_deigenstrain[0])[_qp] = _rotated_cubic_parent_electrostrictive_strain;
  }
  if (_deigenstrain[1])
  {
    _vals[0] = (1.0/2.0)*(_Q11[_qp]*_polar_x[_qp]*std::sin(4*_twist[_qp]) - _Q11[_qp]*_polar_y[_qp]*std::cos(4*_twist[_qp]) + _Q11[_qp]*_polar_y[_qp] - _Q12[_qp]*_polar_x[_qp]*std::sin(4*_twist[_qp]) + _Q12[_qp]*_polar_y[_qp]*std::cos(4*_twist[_qp]) + 3*_Q12[_qp]*_polar_y[_qp] - _Q44[_qp]*_polar_x[_qp]*std::sin(4*_twist[_qp]) + _Q44[_qp]*_polar_y[_qp]*std::cos(4*_twist[_qp]) - _Q44[_qp]*_polar_y[_qp]);  // d eps_1
    _vals[1] = (1.0/2.0)*(-_Q11[_qp]*_polar_x[_qp]*std::sin(4*_twist[_qp]) + _Q11[_qp]*_polar_y[_qp]*std::cos(4*_twist[_qp]) + 3*_Q11[_qp]*_polar_y[_qp] + _Q12[_qp]*_polar_x[_qp]*std::sin(4*_twist[_qp]) - _Q12[_qp]*_polar_y[_qp]*std::cos(4*_twist[_qp]) + _Q12[_qp]*_polar_y[_qp] + _Q44[_qp]*_polar_x[_qp]*std::sin(4*_twist[_qp]) - _Q44[_qp]*_polar_y[_qp]*std::cos(4*_twist[_qp]) + _Q44[_qp]*_polar_y[_qp]);  // d eps_2
    _vals[2] = 2*_Q12[_qp]*_polar_y[_qp];  // d eps_3
    _vals[3] = _Q44[_qp]*_polar_z[_qp];  // d eps_4   23
    _vals[4] = 0;  // d eps_5   13
    _vals[5] = (1.0/2.0)*(-_Q11[_qp]*_polar_x[_qp]*std::cos(4*_twist[_qp]) + _Q11[_qp]*_polar_x[_qp] - _Q11[_qp]*_polar_y[_qp]*std::sin(4*_twist[_qp]) + _Q12[_qp]*_polar_x[_qp]*std::cos(4*_twist[_qp]) - _Q12[_qp]*_polar_x[_qp] + _Q12[_qp]*_polar_y[_qp]*std::sin(4*_twist[_qp]) + _Q44[_qp]*_polar_x[_qp]*std::cos(4*_twist[_qp]) + _Q44[_qp]*_polar_x[_qp] + _Q44[_qp]*_polar_y[_qp]*std::sin(4*_twist[_qp]));  // d eps_6   12
    _rotated_cubic_parent_electrostrictive_strain.fillFromInputVector(_vals);
    (*_deigenstrain[1])[_qp] = _rotated_cubic_parent_electrostrictive_strain;
  }
  if (_deigenstrain[2])
  {
    _vals[0] = 2*_Q12[_qp]*_polar_z[_qp];  // d eps_1
    _vals[1] = 2*_Q12[_qp]*_polar_z[_qp];  // d eps_2
    _vals[2] = 2*_Q11[_qp]*_polar_z[_qp];  // d eps_3
    _vals[3] = _Q44[_qp]*_polar_y[_qp];  // d eps_4   23
    _vals[4] = _Q44[_qp]*_polar_x[_qp];  // d eps_5   13
    _vals[5] = 0;  // d eps_6   12
    _rotated_cubic_parent_electrostrictive_strain.fillFromInputVector(_vals);
    (*_deigenstrain[2])[_qp] = _rotated_cubic_parent_electrostrictive_strain;
  }
}
