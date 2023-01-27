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

#include "FilmSurfaceStressBC.h"

registerMooseObject("FerretApp", FilmSurfaceStressBC);

InputParameters FilmSurfaceStressBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription(
      "stress free surface condition, testing, only for z direction");
  params.addRequiredParam<int>("component","Which component(0 for x, 1 for y, 2 for z) in traction is used");
  params.addRequiredCoupledVar("u_x", "The x component of the local displacement");
  params.addRequiredCoupledVar("u_y", "The y component of the local displacement");
  params.addCoupledVar("u_z", 0.0, "The z component of the local displacement");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

FilmSurfaceStressBC::FilmSurfaceStressBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _component(getParam<int>("component")),
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
    _C11(getMaterialProperty<Real>("C11")),
    _C12(getMaterialProperty<Real>("C12")),
    _C44(getMaterialProperty<Real>("C44")),
    _Q11(getMaterialProperty<Real>("Q11")),
    _Q12(getMaterialProperty<Real>("Q12")),
    _Q44(getMaterialProperty<Real>("Q44"))
{
}

//this is a rather weak boundary condition

Real
FilmSurfaceStressBC::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (_C44[_qp]*(-_polar_x[_qp]*_polar_z[_qp]*_Q44[_qp] + _u_x_grad[_qp](2)) + _C44[_qp]*(-_polar_x[_qp]*_polar_z[_qp]*_Q44[_qp] + _u_z_grad[_qp](0)));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (_C44[_qp]*(-_polar_y[_qp]*_polar_z[_qp]*_Q44[_qp] + _u_y_grad[_qp](2)) + _C44[_qp]*(-_polar_y[_qp]*_polar_z[_qp]*_Q44[_qp] + _u_z_grad[_qp](1)));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (_C12[_qp]*(-_polar_x[_qp]*_polar_x[_qp]*_Q11[_qp] - (_polar_y[_qp]*_polar_y[_qp] + _polar_z[_qp]*_polar_z[_qp])*_Q12[_qp] + _u_x_grad[_qp](0)) + 
 _C12[_qp]*(-_polar_y[_qp]*_polar_y[_qp]*_Q11[_qp] - (_polar_x[_qp]*_polar_x[_qp] + _polar_z[_qp]*_polar_z[_qp])*_Q12[_qp] + _u_z_grad[_qp](2)) + _C11[_qp]*(-_polar_z[_qp]*_polar_z[_qp]*_Q11[_qp] - (_polar_x[_qp]*_polar_x[_qp] + _polar_y[_qp]*_polar_y[_qp])*_Q12[_qp] + _u_z_grad[_qp](2)));
  }
  else
    return 0.0;
}
  
Real
FilmSurfaceStressBC::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _grad_phi[_j][_qp](2) * _C44[_qp];
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _grad_phi[_j][_qp](2) * _C44[_qp];
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _grad_phi[_j][_qp](2) * _C11[_qp];
  }
  else
    return 0.0;
}
