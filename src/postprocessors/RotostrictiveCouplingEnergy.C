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

#include "RotostrictiveCouplingEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", RotostrictiveCouplingEnergy);

InputParameters RotostrictiveCouplingEnergy::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over the coupling free energy density between the AFD and elastic fields.");
  params.addRequiredCoupledVar("u_x", "The x component of the local elastic displacement");
  params.addRequiredCoupledVar("u_y", "The y component of the local elastic displacement");
  params.addCoupledVar("u_z", 0.0, "The z component of the local elastic displacement");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase tilt vector");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase tilt vector");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the antiphase tilt vector");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

RotostrictiveCouplingEnergy::RotostrictiveCouplingEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _u_x_var(coupled("u_x")),
   _u_y_var(coupled("u_y")),
   _u_z_var(coupled("u_z")),
   _antiphase_A_x_var(coupled("antiphase_A_x")),
   _antiphase_A_y_var(coupled("antiphase_A_y")),
   _antiphase_A_z_var(coupled("antiphase_A_z")),
   _u_x_grad(coupledGradient("u_x")),
   _u_y_grad(coupledGradient("u_y")),
   _u_z_grad(coupledGradient("u_z")),
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
   _r11(getMaterialProperty<Real>("r11")),
   _r12(getMaterialProperty<Real>("r12")),
   _r44(getMaterialProperty<Real>("r44")),
   _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
RotostrictiveCouplingEnergy::computeQpIntegral()
{
  return _energy_scale*(-(-2.0*_r44[_qp]*((_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0)))/2.0 + (_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*(_u_x_grad[_qp](2) + _u_z_grad[_qp](0)))/2.0 + (_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)))/2.0) - _r12[_qp]*((Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_u_x_grad[_qp](0) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_u_y_grad[_qp](1) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*_u_z_grad[_qp](2)) - 
   _r11[_qp]*(Utility::pow<2>(_antiphase_A_x[_qp])*_u_x_grad[_qp](0) + Utility::pow<2>(_antiphase_A_y[_qp])*_u_y_grad[_qp](1) + Utility::pow<2>(_antiphase_A_z[_qp])*_u_z_grad[_qp](2))));

}
