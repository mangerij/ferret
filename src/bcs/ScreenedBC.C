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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "ScreenedBC.h"

template<>
InputParameters validParams<ScreenedBC>()
{
    InputParameters params = validParams<IntegratedBC>();
    params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
    params.addRequiredParam<Real>("permittivity", "permittivity");
    params.addRequiredParam<Real>("lambda", "lambda");
    params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
    params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
    params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return params;
}

ScreenedBC::ScreenedBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _potential_int_grad(coupledGradient("potential_int")),
  _permittivity(getParam<Real>("permittivity")),
  _lambda(getParam<Real>("lambda")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z"))
{
}

Real
ScreenedBC::computeQpResidual()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);

  return _normals[_qp] * (_permittivity * _potential_int_grad[_qp] + _lambda * w);
}
