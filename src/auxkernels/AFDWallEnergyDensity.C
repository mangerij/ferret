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

#include "AFDWallEnergyDensity.h"

template<>
InputParameters validParams<AFDWallEnergyDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculates the free energy density due to the local gradients in the antiferrodistortive vector field");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the afd vector field");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the afd vector field");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the afd vector field");
  params.addRequiredParam<Real>("H110","antiphase penalty coefficients");
  params.addRequiredParam<Real>("H11_H110","Ratio of antiphase penalty coefficients");
  params.addRequiredParam<Real>("H12_H110","Ratio of antiphase penalty coefficients");
  params.addRequiredParam<Real>("H44_H110","Ratio of antiphase penalty coefficients");
  params.addRequiredParam<Real>("H44P_H110","Ratio of antiphase penalty coefficients");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

AFDWallEnergyDensity::AFDWallEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
  _antiferrodis_A_x_grad(coupledGradient("antiferrodis_A_x")),
  _antiferrodis_A_y_grad(coupledGradient("antiferrodis_A_y")),
  _antiferrodis_A_z_grad(coupledGradient("antiferrodis_A_z")),
  _H110(getParam<Real>("H110")),
  _H11(getParam<Real>("H11_H110")*_H110),
  _H12(getParam<Real>("H12_H110")*_H110),
  _H44(getParam<Real>("H44_H110")*_H110),
  _H44P(getParam<Real>("H44P_H110")*_H110),
  _len_scale(getParam<Real>("len_scale"))
{}

Real
AFDWallEnergyDensity::computeValue()
{
  return (0.5*_H11*(pow(_antiferrodis_A_x_grad[_qp](0),2)+pow(_antiferrodis_A_y_grad[_qp](1),2.0)+pow(_antiferrodis_A_z_grad[_qp](2),2.0))+
    _H12*(_antiferrodis_A_x_grad[_qp](0)*_antiferrodis_A_y_grad[_qp](1)+_antiferrodis_A_y_grad[_qp](1)*_antiferrodis_A_z_grad[_qp](2)+_antiferrodis_A_x_grad[_qp](0)*_antiferrodis_A_z_grad[_qp](2))+
    0.5*_H44*(pow(_antiferrodis_A_x_grad[_qp](1)+_antiferrodis_A_y_grad[_qp](0),2.0)+pow(_antiferrodis_A_y_grad[_qp](2)+_antiferrodis_A_z_grad[_qp](1),2.0)+pow(_antiferrodis_A_x_grad[_qp](2)+_antiferrodis_A_z_grad[_qp](0),2.0))+ 0.5*_H44P*(pow(_antiferrodis_A_x_grad[_qp](1)-_antiferrodis_A_y_grad[_qp](0),2)+pow(_antiferrodis_A_y_grad[_qp](2)-_antiferrodis_A_z_grad[_qp](1),2.0)+pow(_antiferrodis_A_x_grad[_qp](2)-_antiferrodis_A_z_grad[_qp](0),2)))*_len_scale;
}
