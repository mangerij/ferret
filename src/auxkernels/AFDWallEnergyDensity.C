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

#include "AFDWallEnergyDensity.h"

registerMooseObject("FerretApp", AFDWallEnergyDensity);

InputParameters AFDWallEnergyDensity::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the free energy density due to the local gradients in the antiphasetortive vector field");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase tilt vector field");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase tilt vector field");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the antiphase tilt vector field");
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
  _antiphase_A_x_grad(coupledGradient("antiphase_A_x")),
  _antiphase_A_y_grad(coupledGradient("antiphase_A_y")),
  _antiphase_A_z_grad(coupledGradient("antiphase_A_z")),
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
  return (0.5*_H11*(pow(_antiphase_A_x_grad[_qp](0),2)+pow(_antiphase_A_y_grad[_qp](1),2.0)+pow(_antiphase_A_z_grad[_qp](2),2.0))+
    _H12*(_antiphase_A_x_grad[_qp](0)*_antiphase_A_y_grad[_qp](1)+_antiphase_A_y_grad[_qp](1)*_antiphase_A_z_grad[_qp](2)+_antiphase_A_x_grad[_qp](0)*_antiphase_A_z_grad[_qp](2))+
    0.5*_H44*(pow(_antiphase_A_x_grad[_qp](1)+_antiphase_A_y_grad[_qp](0),2.0)+pow(_antiphase_A_y_grad[_qp](2)+_antiphase_A_z_grad[_qp](1),2.0)+pow(_antiphase_A_x_grad[_qp](2)+_antiphase_A_z_grad[_qp](0),2.0))+ 0.5*_H44P*(pow(_antiphase_A_x_grad[_qp](1)-_antiphase_A_y_grad[_qp](0),2)+pow(_antiphase_A_y_grad[_qp](2)-_antiphase_A_z_grad[_qp](1),2.0)+pow(_antiphase_A_x_grad[_qp](2)-_antiphase_A_z_grad[_qp](0),2)))*_len_scale;
}
