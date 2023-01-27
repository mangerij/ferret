/*
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the ter___Ms of the GNU General Public License as published by
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
#include "LocalConservedLangevinNoise.h"

registerMooseObject("FerretApp", LocalConservedLangevinNoise);

InputParameters
LocalConservedLangevinNoise::validParams()
{
  InputParameters params = LocalLangevinNoise::validParams();
  params.addClassDescription("Source term for noise from a ConservedNoise userobject");
  params.addRequiredParam<UserObjectName>(
      "noise", "ConservedNoise userobject that produces the random numbers");
  params.addRequiredCoupledVar("T_bath", "The localized thermal bath");
  return params;
}
LocalConservedLangevinNoise::LocalConservedLangevinNoise(const InputParameters & parameters)
  : LocalLangevinNoise(parameters),
 _noise(getUserObject<LocalConservedNoiseInterface>("noise")),
  _T_bath(coupledValue("T_bath"))
{
  if (parameters.isParamSetByUser("seed"))
    paramError(
        "seed",
        "This parameter has no effect in this kernel. The noise is generated in the user object "
        "specified in the 'noise' parameter. Specify a seed in that user object instead.");
}

Real
LocalConservedLangevinNoise::computeQpResidual()
{
  return -_T_bath[_qp]*_test[_i][_qp] * _noise.getQpValue(_current_elem->id(), _qp) * _amplitude *
         _multiplier_prop[_qp];
}
