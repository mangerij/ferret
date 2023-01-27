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

//* This file is modified from part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LocalLangevinNoise.h"
#include "MooseRandom.h"

registerMooseObject("FerretApp", LocalLangevinNoise);

InputParameters
LocalLangevinNoise::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Source term for non-conserved Langevin noise");
  params.addRequiredParam<Real>("amplitude", "Amplitude"); // per sqrt(time)");
  params.addParam<MaterialPropertyName>(
      "multiplier",
      1.0,
      "Material property to multiply the random numbers with (defaults to 1.0 if omitted)");
  return params;
}
LocalLangevinNoise::LocalLangevinNoise(const InputParameters & parameters)
  : Kernel(parameters),
    _amplitude(getParam<Real>("amplitude")),
    _multiplier_prop(getMaterialProperty<Real>("multiplier"))
{
}

void
LocalLangevinNoise::residualSetup()
{
  unsigned int rseed = _t_step;
  MooseRandom::seed(rseed);
}

Real
LocalLangevinNoise::computeQpResidual()
{
  return -_test[_i][_qp] * (2.0 * MooseRandom::rand() - 1.0) * _amplitude * _multiplier_prop[_qp];
}
