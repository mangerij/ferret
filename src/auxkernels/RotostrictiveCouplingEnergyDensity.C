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

#include "RotostrictiveCouplingEnergyDensity.h"

registerMooseObject("FerretApp", RotostrictiveCouplingEnergyDensity);

InputParameters RotostrictiveCouplingEnergyDensity::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("disp_x", "The x component of the displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the displacement");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase tilt vector");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase tilt vector");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the antiphase tilt vector");
  params.addRequiredParam<Real>("r11", "The coupling constants");
  params.addRequiredParam<Real>("r12", "The coupling constants");
  params.addRequiredParam<Real>("r44", "The coupling constants");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

RotostrictiveCouplingEnergyDensity::RotostrictiveCouplingEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
   _disp_x_var(coupled("disp_x")),
   _disp_y_var(coupled("disp_y")),
   _disp_z_var(coupled("disp_z")),
   _antiphase_A_x_var(coupled("antiphase_A_x")),
   _antiphase_A_y_var(coupled("antiphase_A_y")),
   _antiphase_A_z_var(coupled("antiphase_A_z")),
   _disp_x_grad(coupledGradient("disp_x")),
   _disp_y_grad(coupledGradient("disp_y")),
   _disp_z_grad(coupledGradient("disp_z")),
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
   _r11(getParam<Real>("r11")),
   _r12(getParam<Real>("r12")),
   _r44(getParam<Real>("r44")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
RotostrictiveCouplingEnergyDensity::computeValue()
{
  return -(-2.0*_r44*((_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0)))/2.0 + (_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*(_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0)))/2.0 + (_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*(_disp_y_grad[_qp](2) + _disp_z_grad[_qp](1)))/2.0) - _r12*((Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_disp_x_grad[_qp](0) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_disp_y_grad[_qp](1) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*_disp_z_grad[_qp](2)) -
   _r11*(Utility::pow<2>(_antiphase_A_x[_qp])*_disp_x_grad[_qp](0) + Utility::pow<2>(_antiphase_A_y[_qp])*_disp_y_grad[_qp](1) + Utility::pow<2>(_antiphase_A_z[_qp])*_disp_z_grad[_qp](2))) * std::pow(_len_scale,3);
}
