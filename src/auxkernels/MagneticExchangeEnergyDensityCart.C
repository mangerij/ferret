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

#include "MagneticExchangeEnergyDensityCart.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MagneticExchangeEnergyDensityCart);

InputParameters MagneticExchangeEnergyDensityCart::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the magnetization");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and aJ");
  return params;
}

MagneticExchangeEnergyDensityCart::MagneticExchangeEnergyDensityCart(const InputParameters & parameters) :
  AuxKernel(parameters),
  _mag_x_grad(coupledGradient("mag_x")),
  _mag_y_grad(coupledGradient("mag_y")),
  _mag_z_grad(coupledGradient("mag_z")),
  _Ae(getMaterialProperty<Real>("Ae")),
  _Ms(getMaterialProperty<Real>("Ms")),
  _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
MagneticExchangeEnergyDensityCart::computeValue()
{
  return _energy_scale*_Ms[_qp]*_Ms[_qp]*(((_Ae[_qp])*(Utility::pow<2>(_mag_x_grad[_qp](0))+Utility::pow<2>(_mag_x_grad[_qp](1))+Utility::pow<2>(_mag_x_grad[_qp](2))+Utility::pow<2>(_mag_y_grad[_qp](0))+Utility::pow<2>(_mag_y_grad[_qp](1))+Utility::pow<2>(_mag_y_grad[_qp](2))+Utility::pow<2>(_mag_z_grad[_qp](0))+Utility::pow<2>(_mag_z_grad[_qp](1))+Utility::pow<2>(_mag_z_grad[_qp](2)))));
}
