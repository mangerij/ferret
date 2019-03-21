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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/


#include "SecondDerivativeImplicitEulerScaled.h"
#include "SubProblem.h"

template <>
InputParameters
validParams<SecondDerivativeImplicitEulerScaled>()
{
  InputParameters params = validParams<TimeKernel>();
  params.addParam<Real>("dampening", 1.0, "The second order dampening");
  return params;
}

SecondDerivativeImplicitEulerScaled::SecondDerivativeImplicitEulerScaled(const InputParameters & parameters)
  :TimeKernel(parameters),
  _u_old(valueOld()),
  _u_older(valueOlder()),
  _damnpening(getParam<Real>("dampening"))
{
}

Real
SecondDerivativeImplicitEulerScaled::computeQpResidual()
{
  return _test[_i][_qp] * _damnpening *((_u[_qp] - 2 * _u_old[_qp] + _u_older[_qp]) / (_dt * _dt));
}

Real
SecondDerivativeImplicitEulerScaled::computeQpJacobian()
{
  return _test[_i][_qp] * _damnpening* (_phi[_j][_qp] / (_dt * _dt));
}
