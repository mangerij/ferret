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

#include "AFMExchangeStiffnessEnergy.h"

registerMooseObject("FerretApp", AFMExchangeStiffnessEnergy);

InputParameters AFMExchangeStiffnessEnergy::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over the DM interaction free energy density (coupling AFD and magnetic ordering).");
  params.addRequiredCoupledVar("Neel_L_x", "The x component of the AFM Neel vector");
  params.addRequiredCoupledVar("Neel_L_y", "The y component of the AFM Neel vector");
  params.addRequiredCoupledVar("Neel_L_z", "The z component of the AFM Neel vector");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV to/from aJ");
  return params;
}

AFMExchangeStiffnessEnergy::AFMExchangeStiffnessEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _Neel_L_x_grad(coupledGradient("Neel_L_x")),
   _Neel_L_y_grad(coupledGradient("Neel_L_y")),
   _Neel_L_z_grad(coupledGradient("Neel_L_z")),
   _Ae(getMaterialProperty<Real>("Ae")),
   _Ms(getMaterialProperty<Real>("Ms")),
   _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
AFMExchangeStiffnessEnergy::computeQpIntegral()
{

//Note that this is natively given in 1/1000 aJ since Neel L is normalized. To convert to eV, we use energy_scale = 0.0187245.

  return _energy_scale*(_Ae[_qp]*(_Neel_L_x_grad[_qp](0)*_Neel_L_x_grad[_qp](0)+_Neel_L_x_grad[_qp](1)*_Neel_L_x_grad[_qp](1)+_Neel_L_x_grad[_qp](2)*_Neel_L_x_grad[_qp](2)+_Neel_L_y_grad[_qp](0)*_Neel_L_y_grad[_qp](0)+_Neel_L_y_grad[_qp](1)*_Neel_L_y_grad[_qp](1)+_Neel_L_y_grad[_qp](2)*_Neel_L_y_grad[_qp](2)+_Neel_L_z_grad[_qp](0)*_Neel_L_z_grad[_qp](0)+_Neel_L_z_grad[_qp](1)*_Neel_L_z_grad[_qp](1)+_Neel_L_z_grad[_qp](2)*_Neel_L_z_grad[_qp](2)));
}
