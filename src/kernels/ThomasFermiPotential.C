/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#include "ThomasFermiPotential.h"

class ThomasFermiPotential;

template<>
InputParameters validParams<ThomasFermiPotential>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("potential_int", "The electrostatic potential to be screened");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  params.addParam<Real>("TFconstant", 0.0, "the Thomas-Fermi constant which proportional to rho divided by the Fermi level");
  return params;
}

ThomasFermiPotential::ThomasFermiPotential(const InputParameters & parameters)
  :Kernel(parameters),
   _potential_int_var(coupled("potential_int")),
   _potential_int(coupledValue("potential_int")),
   _len_scale(getParam<Real>("len_scale")),
   _TFconstant(getParam<Real>("TFconstant"))
{
}

Real
ThomasFermiPotential::computeQpResidual()
{
  return _TFconstant * _potential_int[_qp] * _test[_i][_qp];
}
Real
ThomasFermiPotential::computeQpJacobian()
{
  return _TFconstant * _potential_int[_qp] * _phi[_j][_qp] * _test[_i][_qp];
}

