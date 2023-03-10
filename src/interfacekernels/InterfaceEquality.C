//* This file is part of the MOOSE framework
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
