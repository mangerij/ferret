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

template<>
InputParameters validParams<MagneticExchangeEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over the magnetic exchange energy density.");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetization vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetization vector");
  params.addCoupledVar("mag_z", 0.0, "The z component of the magnetization vector");
  params.addRequiredParam<Real>("A", "The constant of magnetic exchange");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

MagneticExchangeEnergy::MagneticExchangeEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _mag_x_grad(coupledGradient("mag_x")),
  _mag_y_grad(coupledGradient("mag_y")),
  _mag_z_grad(coupledGradient("mag_z")),
  _A(getParam<Real>("A")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
MagneticExchangeEnergy::computeQpIntegral()
{
  return (_A*(std::pow(_mag_x_grad[_qp](0),2) + std::pow(_mag_x_grad[_qp](1),2) + std::pow(_mag_x_grad[_qp](2),2) + std::pow(_mag_y_grad[_qp](0),2) + std::pow(_mag_y_grad[_qp](1),2) + std::pow(_mag_y_grad[_qp](2),2) + std::pow(_mag_z_grad[_qp](0),2) + std::pow(_mag_z_grad[_qp](1),2) + std::pow(_mag_z_grad[_qp](2),2))) * std::pow(_len_scale,1.0);
}
