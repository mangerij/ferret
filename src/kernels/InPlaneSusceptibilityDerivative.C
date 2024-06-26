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

#include "InPlaneSusceptibilityDerivative.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", InPlaneSusceptibilityDerivative);

InputParameters InPlaneSusceptibilityDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates the residual for the local free energy which is an eighth order expansion in the polarization.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  return params;
}

InPlaneSusceptibilityDerivative::InPlaneSusceptibilityDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _polar_x(coupledValue("polar_x")),
   _chi(getMaterialProperty<Real>("chi"))
{
}

Real
InPlaneSusceptibilityDerivative::computeQpResidual()
{
  return _test[_i][_qp]*(1/_chi[_qp])*_polar_x[_qp];
}

Real
InPlaneSusceptibilityDerivative::computeQpJacobian()
{
  return _phi[_j][_qp]*(1/_chi[_qp]);
}

