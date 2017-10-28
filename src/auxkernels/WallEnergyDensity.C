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

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "WallEnergyDensity.h"

template<>
InputParameters validParams<WallEnergyDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("G110","Domain wall penalty coefficients");
  params.addRequiredParam<Real>("G11_G110","Ratio of domain wall penalty coefficients");
  params.addRequiredParam<Real>("G12_G110","Ratio of domain wall penalty coefficients");
  params.addRequiredParam<Real>("G44_G110","Ratio of domain wall penalty coefficients");
  params.addRequiredParam<Real>("G44P_G110","Ratio of domain wall penalty coefficients");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

WallEnergyDensity::WallEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z")),
  _G110(getParam<Real>("G110")),
  _G11(getParam<Real>("G11_G110")*_G110),
  _G12(getParam<Real>("G12_G110")*_G110),
  _G44(getParam<Real>("G44_G110")*_G110),
  _G44P(getParam<Real>("G44P_G110")*_G110),
  _len_scale(getParam<Real>("len_scale"))
{}

Real
WallEnergyDensity::computeValue()
{
  return (0.5*_G11*(pow(_polar_x_grad[_qp](0),2)+pow(_polar_y_grad[_qp](1),2)+pow(_polar_z_grad[_qp](2),2))+
    _G12*(_polar_x_grad[_qp](0)*_polar_y_grad[_qp](1)+_polar_y_grad[_qp](1)*_polar_z_grad[_qp](2)+_polar_x_grad[_qp](0)*_polar_z_grad[_qp](2))+
    0.5*_G44*(pow(_polar_x_grad[_qp](1)+_polar_y_grad[_qp](0),2)+pow(_polar_y_grad[_qp](2)+_polar_z_grad[_qp](1),2)+pow(_polar_x_grad[_qp](2)+_polar_z_grad[_qp](0),2))+
    0.5*_G44P*(pow(_polar_x_grad[_qp](1)-_polar_y_grad[_qp](0),2)+pow(_polar_y_grad[_qp](2)-_polar_z_grad[_qp](1),2)+pow(_polar_x_grad[_qp](2)-_polar_z_grad[_qp](0),2)))*_len_scale;
}
