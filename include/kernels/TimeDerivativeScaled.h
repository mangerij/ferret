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

#ifndef TIMEDERIVATIVESCALED_H
#define TIMEDERIVATIVESCALED_H

#include "TimeKernel.h"
#include "libmesh/quadrature.h"
#include "Assembly.h"

class TimeDerivativeScaled : public TimeKernel
{
public:
  TimeDerivativeScaled(const InputParameters & parameters);

  static InputParameters validParams();

  virtual void computeJacobian();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  bool _lumping;
  const Real _time_scale;
};

#endif //TIMEDERIVATIVESCALED_H
