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

#include "UncoupledGrainBoundaryPhaseIndicator.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", UncoupledGrainBoundaryPhaseIndicator);

InputParameters UncoupledGrainBoundaryPhaseIndicator::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates the residual for the local free energy regarding the GB phase indicator variable");
  params.addRequiredCoupledVar("gbpf", "The grain boundary phase indicator");
  return params;
}

UncoupledGrainBoundaryPhaseIndicator::UncoupledGrainBoundaryPhaseIndicator(const InputParameters & parameters)
  :Kernel(parameters),
   _gbpf_var(coupled("gbpf")),
   _gbpf(coupledValue("gbpf")),
   _a(getMaterialProperty<Real>("a")),
   _b(getMaterialProperty<Real>("b"))
{
}

Real
UncoupledGrainBoundaryPhaseIndicator::computeQpResidual()
{
  return _test[_i][_qp] * (2.0*_a[_qp]*_gbpf[_qp] + 4.0*_b[_qp]*_gbpf[_qp]*_gbpf[_qp]*_gbpf[_qp]);
}

Real
UncoupledGrainBoundaryPhaseIndicator::computeQpJacobian()
{
  return _test[_i][_qp] * _phi[_j][_qp] * (2.0*_a[_qp] + 12.0*_b[_qp]*_gbpf[_qp]*_gbpf[_qp]);
}
