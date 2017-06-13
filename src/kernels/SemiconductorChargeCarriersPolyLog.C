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

#include "SemiconductorChargeCarriersPolyLog.h"

class SemiconductorChargeCarriersPolyLog;

template<>
InputParameters validParams<SemiconductorChargeCarriersPolyLog>()
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

SemiconductorChargeCarriersPolyLog::SemiconductorChargeCarriersPolyLog(const InputParameters & parameters)
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
SemiconductorChargeCarriersPolyLog::computeQpResidual()
{
  Real nm = _NC * std::pow(_kT, 1.5) * std::exp((-_EC + _EF + _q * _potential_int[_qp]) / _kT );
  Real pp = _NV * std::pow(_kT, 1.5) * std::exp((_EF - _EV + _q * _potential_int[_qp]) / _kT );
  Real NAm = _NA;
  Real rho = _q * ( - nm + pp - NAm);
  return rho * _test[_i][_qp]; //might be off by a minus sign...
}
Real
SemiconductorChargeCarriersPolyLog::computeQpJacobian()
{
  return _phi[_j][_qp] * _test[_i][_qp] * std::pow(_kT, 0.5) * (_NC * std::exp((-_EC + _EF + _q * _potential_int[_qp]) / _kT ) + _NV * std::exp((_EF - _EV + _q * _potential_int[_qp]) / _kT ) );
}
