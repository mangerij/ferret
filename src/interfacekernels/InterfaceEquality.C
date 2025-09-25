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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

//* This file is derived from a part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceEquality.h"

registerMooseObject("MooseApp", InterfaceEquality);

InputParameters
InterfaceEquality::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addParam<MaterialPropertyName>("permittivity", "permittivity", "The permittivity.");
  params.addParam<MaterialPropertyName>(
      "permittivity_neighbor", "permittivity_neighbor", "The neighboring permittivity.");
  params.addClassDescription(
      "The kernel is utilized to establish equivalence on an interface for variables.");
  return params;
}

InterfaceEquality::InterfaceEquality(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _permittivity(getMaterialProperty<Real>("permittivity")),
    _permittivity_neighbor(getNeighborMaterialProperty<Real>("permittivity_neighbor"))
{
}

Real
InterfaceEquality::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      r = _test[_i][_qp] * -_permittivity_neighbor[_qp] * _neighbor_value[_qp] ;
      break;

    case Moose::Neighbor:
      r = _test_neighbor[_i][_qp] * _permittivity[_qp] * _u[_qp] ;
      break;
  }

  return r;
}

Real
InterfaceEquality::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;

  switch (type)
  {
    case Moose::ElementElement:
    case Moose::NeighborNeighbor:
      break;

    case Moose::NeighborElement:
      jac = _test_neighbor[_i][_qp] * _permittivity[_qp] * _phi[_j][_qp] ;
      break;

    case Moose::ElementNeighbor:
      jac = _test[_i][_qp] * -_permittivity_neighbor[_qp] * _phi_neighbor[_j][_qp];
      break;
  }

  return jac;
}
