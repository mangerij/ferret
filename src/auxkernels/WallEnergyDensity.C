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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "WallEnergyDensity.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", WallEnergyDensity);

template<>
InputParameters validParams<WallEnergyDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

WallEnergyDensity::WallEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z")),
  _G110(getMaterialProperty<Real>("G110")),
  _G11(getMaterialProperty<Real>("G11_G110")),
  _G12(getMaterialProperty<Real>("G12_G110")),
  _G44(getMaterialProperty<Real>("G44_G110")),
  _G44P(getMaterialProperty<Real>("G44P_G110"))
{}

Real
WallEnergyDensity::computeValue()
{
//note that 1/2 * G_11 should be here, but the Kernel isn't consistent. TODO: To prevent regolding (at this date, Nov 2020) I will remove here.
// To remedy this, a factor of 0.5 can be used on input file param G110 and then multiply by the entered input file param G12 by 2.0...
  return _G110[_qp]*((_G11[_qp]*(Utility::pow<2>(_polar_x_grad[_qp](0))+Utility::pow<2>(_polar_y_grad[_qp](1))+Utility::pow<2>(_polar_z_grad[_qp](2)))+
    _G12[_qp]*(_polar_x_grad[_qp](0)*_polar_y_grad[_qp](1)+_polar_y_grad[_qp](1)*_polar_z_grad[_qp](2)+_polar_x_grad[_qp](0)*_polar_z_grad[_qp](2))+
    _G44[_qp]*(pow(_polar_x_grad[_qp](1)+_polar_y_grad[_qp](0),2)+pow(_polar_y_grad[_qp](2)+_polar_z_grad[_qp](1),2)+pow(_polar_x_grad[_qp](2)+_polar_z_grad[_qp](0),2))+
	  _G44P[_qp]*(Utility::pow<2>(_polar_x_grad[_qp](1)-_polar_y_grad[_qp](0))+Utility::pow<2>(_polar_y_grad[_qp](2)-_polar_z_grad[_qp](1))+Utility::pow<2>(_polar_x_grad[_qp](2)-_polar_z_grad[_qp](0)))));
}
