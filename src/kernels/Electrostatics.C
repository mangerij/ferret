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

#include "Electrostatics.h"

registerMooseObject("FerretApp", Electrostatics);

InputParameters Electrostatics::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  return params;
}

Electrostatics::Electrostatics(const InputParameters & parameters)
  :Kernel(parameters),
   _permittivity(getMaterialProperty<Real>("permittivity"))
{
}

Real
Electrostatics::computeQpResidual()
{
  Real Relec = 0.0;
  Relec += _permittivity[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
  ///  Moose::out << "\n R_elec-"; std::cout << " = " << Relec;
  return Relec;
}

Real
Electrostatics::computeQpJacobian()
{
   return _permittivity[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
