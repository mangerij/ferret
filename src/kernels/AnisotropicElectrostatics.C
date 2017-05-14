/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#include "AnisotropicElectrostatics.h"

template<>
InputParameters validParams<AnisotropicElectrostatics>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("inplane_permittivity", "in plane permittivity");
  params.addRequiredParam<Real>("outofplane_permittivity", "out of plane permittivity");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

AnisotropicElectrostatics::AnisotropicElectrostatics(const InputParameters & parameters)
  :Kernel(parameters),
   _inplane_permittivity(getParam<Real>("inplane_permittivity")),
   _outofplane_permittivity(getParam<Real>("outofplane_permittivity")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
AnisotropicElectrostatics::computeQpResidual()
{
  return _inplane_permittivity * (_grad_u[_qp](0) * _grad_test[_i][_qp](0) + (_outofplane_permittivity/_inplane_permittivity) * _grad_u[_qp](1) * _grad_test[_i][_qp](1) + _grad_u[_qp](2) * _grad_test[_i][_qp](2) ) * _len_scale;
}

Real
AnisotropicElectrostatics::computeQpJacobian()
{
  return _inplane_permittivity * (_grad_phi[_j][_qp](0) * _grad_test[_i][_qp](0) + (_outofplane_permittivity/_inplane_permittivity) * _grad_phi[_j][_qp](1) * _grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2) * _grad_test[_i][_qp](2) ) * _len_scale;
}
