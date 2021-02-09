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

#include "InducedPWallEnergyDensity.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", InducedPWallEnergyDensity);

template<>
InputParameters validParams<InducedPWallEnergyDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("induced_polar_x", "The x component of the induced polarization");
  params.addRequiredCoupledVar("induced_polar_y", "The y component of the induced polarization");
  params.addCoupledVar("induced_polar_z", 0.0, "The z component of the induced polarization");
  params.addRequiredParam<Real>("G110","Domain wall penalty coefficients");
  params.addRequiredParam<Real>("G11_G110","Ratio of domain wall penalty coefficients");
  params.addRequiredParam<Real>("G12_G110","Ratio of domain wall penalty coefficients");
  params.addRequiredParam<Real>("G44_G110","Ratio of domain wall penalty coefficients");
  params.addRequiredParam<Real>("G44P_G110","Ratio of domain wall penalty coefficients");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

InducedPWallEnergyDensity::InducedPWallEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z")),
  _induced_polar_x_grad(coupledGradient("induced_polar_x")),
  _induced_polar_y_grad(coupledGradient("induced_polar_y")),
  _induced_polar_z_grad(coupledGradient("induced_polar_z")),
  _G110(getParam<Real>("G110")),
  _G11(getParam<Real>("G11_G110")*_G110),
  _G12(getParam<Real>("G12_G110")*_G110),
  _G44(getParam<Real>("G44_G110")*_G110),
  _G44P(getParam<Real>("G44P_G110")*_G110),
  _len_scale(getParam<Real>("len_scale"))
{}

Real
InducedPWallEnergyDensity::computeValue()
{
//factors of 2 in this file may different from the kernel or wall energy density itself. TODO: It needs to be carefully checked at a later date. 
//For now [Nov 2020], we just want an order of magnitude estimate of this term 
  return _G11*(_polar_x_grad[_qp](0)*_induced_polar_x_grad[_qp](0)+_polar_y_grad[_qp](1)*_induced_polar_y_grad[_qp](1)+_polar_z_grad[_qp](2)*_induced_polar_z_grad[_qp](2))+_G12*(_polar_z_grad[_qp](2)*(_induced_polar_x_grad[_qp](0)+_induced_polar_y_grad[_qp](1))+_polar_y_grad[_qp](1)*(_induced_polar_x_grad[_qp](0)+_induced_polar_z_grad[_qp](2))+_polar_z_grad[_qp](2)*(_induced_polar_y_grad[_qp](1)+_induced_polar_z_grad[_qp](2)))+_G44*((_polar_x_grad[_qp](1)+_polar_y_grad[_qp](0))*(_induced_polar_x_grad[_qp](1)+_induced_polar_y_grad[_qp](0))+(_polar_x_grad[_qp](2)+_polar_z_grad[_qp](0))*(_induced_polar_x_grad[_qp](2)+_induced_polar_z_grad[_qp](0))+(_polar_y_grad[_qp](2)+_polar_z_grad[_qp](1))*(_induced_polar_y_grad[_qp](2)+_induced_polar_z_grad[_qp](1)));
}
