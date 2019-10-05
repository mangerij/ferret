/**
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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "MagneticExchangeEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MagneticExchangeEnergy);

template<>
InputParameters validParams<MagneticExchangeEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over the magnetic exchange energy density.");
  params.addRequiredCoupledVar("magnetic_x", "The x component of the magnetization");
  params.addRequiredCoupledVar("magnetic_y", "The y component of the magnetization");
  params.addCoupledVar("magnetic_z", 0.0, "The z component of the magnetization");
  params.addRequiredParam<Real>("Ae", "Ae");
  params.addRequiredParam<Real>("Ms", "Ms");
  return params;
}

MagneticExchangeEnergy::MagneticExchangeEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _magnetic_x_grad(coupledGradient("magnetic_x")),
  _magnetic_y_grad(coupledGradient("magnetic_y")),
  _magnetic_z_grad(coupledGradient("magnetic_z")),
  _Ae(getParam<Real>("Ae")),
  _Ms(getParam<Real>("Ms"))
{
}

Real
MagneticExchangeEnergy::computeQpIntegral()
{
  return (_Ae/_Ms)*(Utility::pow<2>(_magnetic_x_grad[_qp](0))+Utility::pow<2>(_magnetic_x_grad[_qp](1))+Utility::pow<2>(_magnetic_x_grad[_qp](2))+Utility::pow<2>(_magnetic_y_grad[_qp](0))+Utility::pow<2>(_magnetic_y_grad[_qp](1))+Utility::pow<2>(_magnetic_y_grad[_qp](2))+Utility::pow<2>(_magnetic_z_grad[_qp](0))+Utility::pow<2>(_magnetic_z_grad[_qp](1))+Utility::pow<2>(_magnetic_z_grad[_qp](2)));
}
