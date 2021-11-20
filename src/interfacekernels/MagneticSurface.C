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

#include "MagneticSurface.h"

#include <cmath>

registerMooseObject("FerretApp", MagneticSurface);

InputParameters MagneticSurface::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  return params;
}

MagneticSurface::MagneticSurface(const InputParameters & parameters) :
    InterfaceKernel(parameters),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _mu0(getMaterialProperty<Real>("mu0")),
  _Ms(getMaterialProperty<Real>("Ms"))
{
  if (!parameters.isParamValid("boundary"))
  {
    mooseError("In order to use the MagneticSurface dgkernel, you must specify a boundary where it will live.");
  }
}

Real
MagneticSurface::computeQpResidual(Moose::DGResidualType type)
{
  Real r = ((_grad_u[_qp] * _normals[_qp] - _grad_neighbor_value[_qp] * _normals[_qp]) - _Ms[_qp]*_mu0[_qp]*(_mag_x[_qp]*_normals[_qp](0)+_mag_y[_qp]*_normals[_qp](1)+_mag_z[_qp]*_normals[_qp](2)));

  switch (type)
  {
  case Moose::Element:
    r *= _test[_i][_qp];
    break;

  case Moose::Neighbor:
    r *= -_test_neighbor[_i][_qp]; //sign comes from "opposite" face, not difference in var
    break;
  }

  return r;
}

Real
MagneticSurface::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;

  switch (type)
  {
    case Moose::ElementElement:
      jac += _grad_phi[_j][_qp] * _normals[_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      jac -= -_grad_phi_neighbor[_j][_qp] * _normals[_qp] * _test_neighbor[_i][_qp];  //first minus comes from residual definition of test_neighbor, second comes from physics
      break;

    case Moose::NeighborElement:
      jac -= _grad_phi[_j][_qp] * _normals[_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::ElementNeighbor:
      jac +=  -_grad_phi_neighbor[_j][_qp] * _normals[_qp] * _test[_i][_qp];
      break;
  }

  return jac;
  //missing off-diagonals but let's zero them for now...
  //signs seem a little off on the jac terms (opposite?)
}
