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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "AnisotropicElectrostatics.h"

registerMooseObject("FerretApp", AnisotropicElectrostatics);

template<>
InputParameters validParams<AnisotropicElectrostatics>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  return params;
}

AnisotropicElectrostatics::AnisotropicElectrostatics(const InputParameters & parameters)
  :Kernel(parameters),
   _eps1(getMaterialProperty<Real>("eps1")),
   _eps2(getMaterialProperty<Real>("eps2")),
   _eps3(getMaterialProperty<Real>("eps3"))
{
}

Real
AnisotropicElectrostatics::computeQpResidual()
{
  Real Relec = 0.0;
  Relec += (_eps1[_qp] * _grad_u[_qp](0) * _grad_test[_i][_qp](0)+ _eps2[_qp] * _grad_u[_qp](1) * _grad_test[_i][_qp](1) +_eps3[_qp] * _grad_u[_qp](2) * _grad_test[_i][_qp](2) );
  ///  Moose::out << "\n R_elec-"; std::cout << " = " << Relec;
  return Relec;
}

Real
AnisotropicElectrostatics::computeQpJacobian()
{
   return (_eps1[_qp] * _grad_phi[_j][_qp](0) * _grad_test[_i][_qp](0)+ _eps2[_qp] * _grad_phi[_j][_qp](1) * _grad_test[_i][_qp](1) +_eps3[_qp] * _grad_phi[_j][_qp](2) * _grad_test[_i][_qp](2) );
}
