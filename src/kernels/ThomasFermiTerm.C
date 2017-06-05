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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "ThomasFermiTerm.h"

class ThomasFermiTerm;

template<>
InputParameters validParams<ThomasFermiTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("potential_int", "The electrostatic potential to be screened");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  params.addParam<Real>("q", 1.0, "the electron charge");
  params.addParam<Real>("rho", 0.0, "bulk conduction electron density");
  params.addParam<Real>("EF", 1.0, "Fermi level");
  return params;
}

ThomasFermiTerm::ThomasFermiTerm(const InputParameters & parameters)
  :Kernel(parameters),
   _potential_int_var(coupled("potential_int")),
   _potential_int(coupledValue("potential_int")),
   _len_scale(getParam<Real>("len_scale")),
   _q(getParam<Real>("q")),
   _rho(getParam<Real>("rho")),
   _EF(getParam<Real>("EF"))
{
}

Real
ThomasFermiTerm::computeQpResidual()
{
  return - _rho * ( std::pow(1.0 - (_q / _EF ) * _potential_int[_qp], 1.5) - 1.0) * _test[_i][_qp];
}
Real
ThomasFermiTerm::computeQpJacobian()
{
  return _rho * (_q / _EF ) * 1.5 * std::pow(1.0 - (_q / _EF ) * _potential_int[_qp], 0.5) * _phi[_j][_qp] * _test[_i][_qp];
}
