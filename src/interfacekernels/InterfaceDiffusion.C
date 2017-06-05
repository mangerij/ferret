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

#include "InterfaceDiffusion.h"

#include <cmath>

template<>
InputParameters validParams<InterfaceDiffusion>()
{
  InputParameters params = validParams<InterfaceKernel>();
  params.addParam<Real>("D",1.,"The diffusion coefficient.");
  params.addParam<Real>("D_neighbor",1.,"The neighboring diffusion coefficient.");
  return params;
}

InterfaceDiffusion::InterfaceDiffusion(const InputParameters & parameters) :
    InterfaceKernel(parameters),
    _D(getParam<Real>("D")),
    _D_neighbor(getParam<Real>("D_neighbor"))
{
  if (!parameters.isParamValid("boundary"))
  {
    mooseError("In order to use the InterfaceDiffusion dgkernel, you must specify a boundary where it will live.");
  }
}

Real
InterfaceDiffusion::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0.5 * (-_D * _grad_u[_qp] * _normals[_qp] + -_D_neighbor * _grad_neighbor_value[_qp] * _normals[_qp]);

  switch (type)
  {
  case Moose::Element:
    r *= _test[_i][_qp];
    break;

  case Moose::Neighbor:
    r *= -_test_neighbor[_i][_qp];
    break;
  }

  return r;
}

Real
InterfaceDiffusion::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;

  switch (type)
  {

    case Moose::ElementElement:
      jac -= 0.5 * _D * _grad_phi[_j][_qp] * _normals[_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      jac += 0.5 * _D_neighbor * _grad_phi_neighbor[_j][_qp] * _normals[_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::NeighborElement:
      jac += 0.5 * _D * _grad_phi[_j][_qp] * _normals[_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::ElementNeighbor:
      jac -=  0.5 * _D_neighbor * _grad_phi_neighbor[_j][_qp] * _normals[_qp] * _test[_i][_qp];
      break;
  }

  return jac;
}
