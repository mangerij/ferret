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

#include "FluctuationKernel.h"
#include<cmath>

registerMooseObject("FerretApp", FluctuationKernel);

InputParameters FluctuationKernel::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution introducing fluctuations useful in quasi-static hysteretic switching.");
  params.addRequiredCoupledVar("deltaPi", "The magnitude of the fluctuation across the ith component");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

FluctuationKernel::FluctuationKernel(const InputParameters & parameters)
  :Kernel(parameters),
   _deltaPi(coupledValue("deltaPi")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
FluctuationKernel::computeQpResidual()
{
  return -_deltaPi[_qp] * _test[_i][_qp];
}
