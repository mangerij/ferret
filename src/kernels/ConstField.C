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

#include "ConstField.h"

class ConstField;

registerMooseObject("FerretApp", ConstField);

InputParameters ConstField::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("This is just a test kernel. It is a residual contribution due to a constant electric field term"
                             " along the z-direction of polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("field", 0.0, "the constant field");
  return params;
}

ConstField::ConstField(const InputParameters & parameters)
  :Kernel(parameters),
   _polar_z_var(coupled("polar_z")),
   _polar_z(coupledValue("polar_z")),
   _len_scale(getParam<Real>("len_scale")),
   _field(getParam<Real>("field"))
{
}

Real
ConstField::computeQpResidual()
{
  return _field * _test[_i][_qp] * std::pow(_len_scale, 3.0);
}

Real
ConstField::computeQpJacobian()
{
  return 0.0;
}
