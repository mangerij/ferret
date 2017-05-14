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

#include "SemiconductorChargeCarriers.h"

class SemiconductorChargeCarriers;

template<>
InputParameters validParams<SemiconductorChargeCarriers>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("potential_int", "The electrostatic potential");
  params.addParam<Real>("q", 1.0, "the variable q corresponding to the value of the charge");
  params.addParam<Real>("kT", 1.0, "kT");
  params.addParam<Real>("NA", 1.0, "NA");
  params.addParam<Real>("NC", 1.0, "NC");
  params.addParam<Real>("NV", 1.0, "NV");
  params.addParam<Real>("EA", 1.0, "EA");
  params.addParam<Real>("EC", 1.0, "EC");
  params.addParam<Real>("EV", 1.0, "EV");
  params.addParam<Real>("EF", 1.0, "EF");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

SemiconductorChargeCarriers::SemiconductorChargeCarriers(const InputParameters & parameters)
  :Kernel(parameters),
   _potential_int_var(coupled("potential_int")),
   _potential_int(coupledValue("potential_int")),
   _q(getParam<Real>("q")),
   _kT(getParam<Real>("kT")),
   _NA(getParam<Real>("NA")),
   _NC(getParam<Real>("NC")),
   _NV(getParam<Real>("NV")),
   _EA(getParam<Real>("EA")),
   _EC(getParam<Real>("EC")),
   _EV(getParam<Real>("EV")),
   _EF(getParam<Real>("EF")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
SemiconductorChargeCarriers::computeQpResidual()
{
  Real nm = _NC * std::pow((std::exp(_q * (_EC - _EF - _potential_int[_qp]) / _kT ) + 1.0), -1.0);
  Real pp = _NV * (1 - std::pow((std::exp(_q * (_EV - _EF - _potential_int[_qp]) / _kT ) + 1.0), -1.0));
  Real NAm = _NA * std::pow((std::exp(_q * (_EA - _EF - _potential_int[_qp]) / _kT ) + 1.0), -1.0);
  Real rho = _q * ( - nm + pp - NAm);
  return rho * _test[_i][_qp]; //might be off by a minus sign...
}
Real
SemiconductorChargeCarriers::computeQpJacobian()
{
  return - _phi[_j][_qp] * _test[_i][_qp] * _q * _q * (1.0 / 4.0*_kT) * (
    _NA * std::pow(1.0/std::cosh(_q * (_EA - _EF - _potential_int[_qp]) / (2.0 * _kT) ), 2.0)
    +_NC * std::pow(1.0/std::cosh(_q * (_EC - _EF - _potential_int[_qp]) / (2.0 * _kT) ), 2.0)
    +_NV * std::pow(1.0/std::cosh(_q * (_EF - _EV + _potential_int[_qp]) / (2.0 * _kT) ), 2.0)
  );
}
