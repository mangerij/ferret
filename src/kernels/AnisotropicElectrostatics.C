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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "AnisotropicElectrostatics.h"

registerMooseObject("FerretApp", AnisotropicElectrostatics);

template<>
InputParameters validParams<AnisotropicElectrostatics>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  params.addRequiredParam<Real>("eps1", "first component of background permitivitty");
  params.addRequiredParam<Real>("eps2", "second component of background permitivitty");
  params.addRequiredParam<Real>("eps3", "third component of background permitivitty");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

AnisotropicElectrostatics::AnisotropicElectrostatics(const InputParameters & parameters)
  :Kernel(parameters),
   _eps1(getParam<Real>("eps1")),
   _eps2(getParam<Real>("eps2")),
   _eps3(getParam<Real>("eps3")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
AnisotropicElectrostatics::computeQpResidual()
{
  Real Relec = 0.0;
  Relec += (_eps1 * _grad_u[_qp](0) * _grad_test[_i][_qp](0)+ _eps2 * _grad_u[_qp](1) * _grad_test[_i][_qp](1) +_eps3 * _grad_u[_qp](2) * _grad_test[_i][_qp](2) )* _len_scale;
  ///  Moose::out << "\n R_elec-"; std::cout << " = " << Relec;
  return Relec;
}

Real
AnisotropicElectrostatics::computeQpJacobian()
{
   return (_eps1 * _grad_phi[_j][_qp](0) * _grad_test[_i][_qp](0)+ _eps2 * _grad_phi[_j][_qp](1) * _grad_test[_i][_qp](1) +_eps3 * _grad_phi[_j][_qp](2) * _grad_test[_i][_qp](2) )* _len_scale;
}
