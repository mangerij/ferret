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

#include "ElectrostrictiveCouplingEnergy.h"

registerMooseObject("FerretApp", ElectrostrictiveCouplingEnergy);

template<>
InputParameters validParams<ElectrostrictiveCouplingEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over the coupling energy density between the elastic and AFD fields.");
  params.addRequiredCoupledVar("u_x", "The x component of the local displacement");
  params.addRequiredCoupledVar("u_y", "The y component of the local displacement");
  params.addCoupledVar("u_z", 0.0, "The z component of the local displacement");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization vector");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization vector");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization vector");
  return params;
}

ElectrostrictiveCouplingEnergy::ElectrostrictiveCouplingEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _u_x_var(coupled("u_x")),
   _u_y_var(coupled("u_y")),
   _u_z_var(coupled("u_z")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _u_x_grad(coupledGradient("u_x")),
   _u_y_grad(coupledGradient("u_y")),
   _u_z_grad(coupledGradient("u_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _q11(getMaterialProperty<Real>("q11")),
   _q12(getMaterialProperty<Real>("q12")),
   _q44(getMaterialProperty<Real>("q44"))
{
}

Real
ElectrostrictiveCouplingEnergy::computeQpIntegral()
{
  return -0.5*(-2.0*_q44[_qp]*((_polar_x[_qp]*_polar_y[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0)))/2.0 + (_polar_x[_qp]*_polar_z[_qp]*(_u_x_grad[_qp](2) + _u_z_grad[_qp](0)))/2.0 + (_polar_y[_qp]*_polar_z[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)))/2.0) - _q12[_qp]*((Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_u_x_grad[_qp](0) + (Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_u_y_grad[_qp](1) + (Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*_u_z_grad[_qp](2)) - 
   _q11[_qp]*(Utility::pow<2>(_polar_x[_qp])*_u_x_grad[_qp](0) + Utility::pow<2>(_polar_y[_qp])*_u_y_grad[_qp](1) + Utility::pow<2>(_polar_z[_qp])*_u_z_grad[_qp](2)));
}
