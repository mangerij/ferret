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

#include "MicroforceElectrostrictiveCouplingEnergy.h"

registerMooseObject("FerretApp", MicroforceElectrostrictiveCouplingEnergy);

template<>
InputParameters validParams<MicroforceElectrostrictiveCouplingEnergy>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Computes the free energy density of the local electrostrictive coupling.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("disp_x", "The x component of the displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the displacement");
  params.addCoupledVar("disp_z", "The z component of the displacement");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization vector");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization vector");
  params.addCoupledVar("polar_z", "The z component of the polarization vector");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

MicroforceElectrostrictiveCouplingEnergy::MicroforceElectrostrictiveCouplingEnergy(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _disp_x_var(coupled("disp_x")),
   _disp_y_var(coupled("disp_y")),
   _disp_z_var(coupled("disp_z")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _disp_x_grad(coupledGradient("disp_x")),
   _disp_y_grad(coupledGradient("disp_y")),
   _disp_z_grad(coupledGradient("disp_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _q11(getMaterialProperty<Real>("q11")),
   _q12(getMaterialProperty<Real>("q12")),
   _q44(getMaterialProperty<Real>("q44")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
MicroforceElectrostrictiveCouplingEnergy::computeValue()
{
  if (_component == 0)
  {
   return -((-2*_polar_x[_qp]*_q11[_qp]*_disp_x_grad[_qp](0) - 2*_q44[_qp]*((_polar_y[_qp]*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0)))/2. + (_polar_z[_qp]*(_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0)))/2.) - _q12[_qp]*(2*_polar_x[_qp]*_disp_y_grad[_qp](1) + 2*_polar_x[_qp]*_disp_z_grad[_qp](2))));
  }
  else if (_component == 1)
  {
    return -((-2*_polar_y[_qp]*_q11[_qp]*_disp_y_grad[_qp](1) - 2*_q44[_qp]*((_polar_x[_qp]*(_disp_x_grad[_qp](1) + _disp_y_grad[_qp](0)))/2. + (_polar_z[_qp]*(_disp_y_grad[_qp](2) + _disp_z_grad[_qp](1)))/2.) - _q12[_qp]*(2*_polar_y[_qp]*_disp_x_grad[_qp](0) + 2*_polar_y[_qp]*_disp_z_grad[_qp](2))));
  }
  else if (_component == 2)
  {
    return -((-(_q12[_qp]*(2*_polar_z[_qp]*_disp_x_grad[_qp](0) + 2*_polar_z[_qp]*_disp_y_grad[_qp](1))) - 2*_q44[_qp]*((_polar_x[_qp]*(_disp_x_grad[_qp](2) + _disp_z_grad[_qp](0)))/2. + (_polar_y[_qp]*(_disp_y_grad[_qp](2) + _disp_z_grad[_qp](1)))/2.) - 2*_polar_z[_qp]*_q11[_qp]*_disp_z_grad[_qp](2)));
  }
  else
    return 0.0;
}