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

#include "RotatedWallEnergy.h"

registerMooseObject("FerretApp", RotatedWallEnergy);

InputParameters RotatedWallEnergy::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over the gradient (Ginzburg) energy of a cubic parent phase, with the crystal axes turned by the material property twist about z (the energy of RotatedWallEnergyDerivative).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

RotatedWallEnergy::RotatedWallEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z")),
  _G11(getMaterialProperty<Real>("G11")),
  _G12(getMaterialProperty<Real>("G12")),
  _G44(getMaterialProperty<Real>("G44")),
  _twist(getMaterialProperty<Real>("twist")),
  _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
RotatedWallEnergy::computeQpIntegral()
{
  return _energy_scale*((3*Utility::pow<2>(_polar_x_grad[_qp](0))*_G11[_qp] + Utility::pow<2>(_polar_x_grad[_qp](1))*_G11[_qp] + 2*_polar_x_grad[_qp](1)*_polar_y_grad[_qp](0)*_G11[_qp] + Utility::pow<2>(_polar_y_grad[_qp](0))*_G11[_qp] + 2*_polar_x_grad[_qp](0)*_polar_y_grad[_qp](1)*_G11[_qp] + 3*Utility::pow<2>(_polar_y_grad[_qp](1))*_G11[_qp] + 4*Utility::pow<2>(_polar_z_grad[_qp](2))*_G11[_qp] + Utility::pow<2>(_polar_x_grad[_qp](0))*_G12[_qp] - Utility::pow<2>(_polar_x_grad[_qp](1))*_G12[_qp] - 
     2*_polar_x_grad[_qp](1)*_polar_y_grad[_qp](0)*_G12[_qp] - Utility::pow<2>(_polar_y_grad[_qp](0))*_G12[_qp] + 6*_polar_x_grad[_qp](0)*_polar_y_grad[_qp](1)*_G12[_qp] + Utility::pow<2>(_polar_y_grad[_qp](1))*_G12[_qp] + 8*_polar_x_grad[_qp](0)*_polar_z_grad[_qp](2)*_G12[_qp] + 8*_polar_y_grad[_qp](1)*_polar_z_grad[_qp](2)*_G12[_qp] + 
     2*(Utility::pow<2>(_polar_x_grad[_qp](1) + _polar_y_grad[_qp](0)) + Utility::pow<2>(_polar_x_grad[_qp](0) - _polar_y_grad[_qp](1)) + 2*(Utility::pow<2>(_polar_x_grad[_qp](2) + _polar_z_grad[_qp](0)) + Utility::pow<2>(_polar_y_grad[_qp](2) + _polar_z_grad[_qp](1))))*_G44[_qp] + 
     (_G11[_qp] - _G12[_qp] - 2*_G44[_qp])*((_polar_x_grad[_qp](0) - _polar_x_grad[_qp](1) - _polar_y_grad[_qp](0) - _polar_y_grad[_qp](1))*(_polar_x_grad[_qp](0) + _polar_x_grad[_qp](1) + _polar_y_grad[_qp](0) - _polar_y_grad[_qp](1))*std::cos(4*_twist[_qp]) + 2*(_polar_x_grad[_qp](1) + _polar_y_grad[_qp](0))*(_polar_x_grad[_qp](0) - _polar_y_grad[_qp](1))*std::sin(4*_twist[_qp])))/8.);
}
