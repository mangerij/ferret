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

#include "DepolEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", DepolEnergy);

InputParameters DepolEnergy::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to an arbitrary depolarization energy adjustment to the P$*$E term.");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<PostprocessorName>("avePz", "the average polarization due to the depol field term");
  params.addParam<Real>("lambda", 1.0, "the screening length term");
  params.addParam<Real>("permitivitty", 1.0, "the permitivitty term");
  return params;
}

DepolEnergy::DepolEnergy(const InputParameters & parameters)
  :Kernel(parameters),
   _polar_z_var(coupled("polar_z")),
   _polar_z(coupledValue("polar_z")),
   _len_scale(getParam<Real>("len_scale")),
   _avePz(getPostprocessorValue("avePz")),
   _lambda(getParam<Real>("lambda")),
   _permitivitty(getParam<Real>("permitivitty"))
{
}

Real
DepolEnergy::computeQpResidual()
{
  return 0.5 * _lambda * (1.0 / _permitivitty) * _avePz * _test[_i][_qp] * Utility::pow<3>(_len_scale);
}

Real
DepolEnergy::computeQpJacobian()
{
  return 0.0;
}
