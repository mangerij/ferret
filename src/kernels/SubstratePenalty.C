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

#include "SubstratePenalty.h"
#include "libmesh/utility.h"

class SubstratePenalty;

registerMooseObject("FerretApp", SubstratePenalty);

InputParameters SubstratePenalty::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates an energy penalty due to deviating from elastic condition imposed by the substrate");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("u_x", "The x component of the local elastic displacement");
  params.addRequiredCoupledVar("u_y", "The y component of the local elastic displacement");
  params.addCoupledVar("u_z", 0.0, "The z component of the local elastic displacement");
  params.addRequiredParam<PostprocessorName>("aveexx", "aveexx");
  params.addRequiredParam<PostprocessorName>("aveeyy", "aveeyy");
  return params;
}

SubstratePenalty::SubstratePenalty(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _aveexx(getPostprocessorValue("aveexx")),
   _aveeyy(getPostprocessorValue("aveeyy")),
   _u_x_grad(coupledGradient("u_x")),
   _u_y_grad(coupledGradient("u_y")),
   _u_z_grad(coupledGradient("u_z")),
   _K(getMaterialProperty<Real>("K")),
   _e0xx(getMaterialProperty<Real>("e0xx")),
   _e0yy(getMaterialProperty<Real>("e0yy"))
{
}

Real
SubstratePenalty::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp]*_K[_qp]*(_u_x_grad[_qp](0) + _aveexx - _e0xx[_qp])*(_u_x_grad[_qp](0) + _aveexx - _e0xx[_qp]);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp]*_K[_qp]*(_u_y_grad[_qp](1) + _aveeyy - _e0yy[_qp])*(_u_y_grad[_qp](1) + _aveeyy - _e0yy[_qp]);
  }
  else
    return 0.0;
}

Real
SubstratePenalty::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp]*_K[_qp]*2.0*_phi[_j][_qp]*(_u_x_grad[_qp](0) + _aveexx - _e0xx[_qp]);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp]*_K[_qp]*2.0*_phi[_j][_qp]*(_u_y_grad[_qp](1) + _aveeyy - _e0yy[_qp]);
  }
  else
    return 0.0;
}
