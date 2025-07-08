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

   You should have received a co_polar_y[_qp] of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ElectrostrictiveCouplingPEnergy.h"

registerMooseObject("FerretApp", ElectrostrictiveCouplingPEnergy);

InputParameters ElectrostrictiveCouplingPEnergy::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates a volume integral over the electrostrictive coupling energy density.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization vector");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization vector");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization vector");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

ElectrostrictiveCouplingPEnergy::ElectrostrictiveCouplingPEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _Q11(getMaterialProperty<Real>("Q11")),
   _Q12(getMaterialProperty<Real>("Q12")),
   _Q44(getMaterialProperty<Real>("Q44")),
   _C11(getMaterialProperty<Real>("C11")),
   _C12(getMaterialProperty<Real>("C12")),
   _C44(getMaterialProperty<Real>("C44")),
   _energy_scale(getParam<Real>("energy_scale")),
   _strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_strain"))
{
}

Real
ElectrostrictiveCouplingPEnergy::computeQpIntegral()
{
  return _energy_scale*(2*_C44[_qp]*(Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*Utility::pow<2>(_Q44[_qp]) + Utility::pow<2>(_strain[_qp](0,1)) + Utility::pow<2>(_strain[_qp](0,2)) - 2*_polar_x[_qp]*_Q44[_qp]*(_polar_y[_qp]*_strain[_qp](0,1) + _polar_z[_qp]*_strain[_qp](0,2)) +
      Utility::pow<2>(-(_polar_y[_qp]*_polar_z[_qp]*_Q44[_qp]) + _strain[_qp](1,2))) + (_C11[_qp]*(Utility::pow<4>(_polar_z[_qp])*Utility::pow<2>(_Q11[_qp]) + 2*Utility::pow<4>(_polar_z[_qp])*Utility::pow<2>(_Q12[_qp]) +
        Utility::pow<4>(_polar_x[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*Utility::pow<2>(_Q12[_qp])) + Utility::pow<4>(_polar_y[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*Utility::pow<2>(_Q12[_qp])) - 2*Utility::pow<2>(_polar_z[_qp])*_Q12[_qp]*_strain[_qp](0,0) +
        Utility::pow<2>(_strain[_qp](0,0)) - 2*Utility::pow<2>(_polar_z[_qp])*_Q12[_qp]*_strain[_qp](1,1) + Utility::pow<2>(_strain[_qp](1,1)) - 2*Utility::pow<2>(_polar_z[_qp])*_Q11[_qp]*_strain[_qp](2,2) + Utility::pow<2>(_strain[_qp](2,2)) +
        2*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_z[_qp])*_Q12[_qp]*(2*_Q11[_qp] + _Q12[_qp]) - _Q11[_qp]*_strain[_qp](1,1) - _Q12[_qp]*(_strain[_qp](0,0) + _strain[_qp](2,2))) +
        2*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp])*_Q12[_qp]*(2*_Q11[_qp] + _Q12[_qp]) + Utility::pow<2>(_polar_z[_qp])*_Q12[_qp]*(2*_Q11[_qp] + _Q12[_qp]) - _Q11[_qp]*_strain[_qp](0,0) - _Q12[_qp]*(_strain[_qp](1,1) + _strain[_qp](2,2)))))/2. +
   _C12[_qp]*(2*Utility::pow<4>(_polar_z[_qp])*_Q11[_qp]*_Q12[_qp] + Utility::pow<4>(_polar_z[_qp])*Utility::pow<2>(_Q12[_qp]) + Utility::pow<4>(_polar_x[_qp])*_Q12[_qp]*(2*_Q11[_qp] + _Q12[_qp]) + Utility::pow<4>(_polar_y[_qp])*_Q12[_qp]*(2*_Q11[_qp] + _Q12[_qp]) -
      Utility::pow<2>(_polar_z[_qp])*_Q11[_qp]*_strain[_qp](0,0) - Utility::pow<2>(_polar_z[_qp])*_Q12[_qp]*_strain[_qp](0,0) - Utility::pow<2>(_polar_z[_qp])*_Q11[_qp]*_strain[_qp](1,1) - Utility::pow<2>(_polar_z[_qp])*_Q12[_qp]*_strain[_qp](1,1) + _strain[_qp](0,0)*_strain[_qp](1,1) -
      2*Utility::pow<2>(_polar_z[_qp])*_Q12[_qp]*_strain[_qp](2,2) + _strain[_qp](0,0)*_strain[_qp](2,2) + _strain[_qp](1,1)*_strain[_qp](2,2) +
      Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp])) + Utility::pow<2>(_polar_z[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp])) -
         _Q11[_qp]*(_strain[_qp](1,1) + _strain[_qp](2,2)) - _Q12[_qp]*(2*_strain[_qp](0,0) + _strain[_qp](1,1) + _strain[_qp](2,2))) +
      Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_z[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp])) - _Q11[_qp]*(_strain[_qp](0,0) + _strain[_qp](2,2)) - _Q12[_qp]*(_strain[_qp](0,0) + 2*_strain[_qp](1,1) + _strain[_qp](2,2)))));
}
