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

#pragma once

#include "LocalConservedNoiseBase.h"
#include "LocalConservedUniformNoiseVeneer.h"

// Forward delcarations

/**
 * Userobject that generates a uniformly distributed random number in the interval [-1:1]
 * once per timestep for every quadrature point in a way that the integral
 * over all random numbers is zero.
 *
 * \see ConservedNoiseBase
 */
class LocalConservedUniformNoise : public LocalConservedUniformNoiseVeneer<LocalConservedNoiseBase>
{
public:
  static InputParameters validParams();

  LocalConservedUniformNoise(const InputParameters & parameters)
    : LocalConservedUniformNoiseVeneer<LocalConservedNoiseBase>(parameters)
  {
  }
};
