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

#include "CarrierRec.h"

registerMooseObject("FerretApp", CarrierRec);

InputParameters
CarrierRec::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to free carriers");
  params.addRequiredParam<Real>("Ev",  "Property name of the Valence band energy (J)");
  params.addRequiredParam<Real>("Ec", "Property name of the Conduction band energy(J)");
  params.addRequiredParam<Real>("Nv","Effective DOS of the valence band(T=298)");
  params.addRequiredParam<Real>("Nc","Effective DOS of the conduction band(T=298)");
  params.addRequiredParam<Real>("T", "temperature (K)");
  params.addRequiredParam<Real>("Kb", "Boltzmann Constant (aJ/K)");
  params.addRequiredParam<Real>("q", "eV (aJ)");
  params.addRequiredParam<Real>("b", "lifetime of carriers");
  return params;
}

CarrierRec::CarrierRec(const InputParameters & parameters)
  :Kernel(parameters),
  _Ev(getParam<Real>("Ev")),
  _Ec(getParam<Real>("Ec")),
  _Nv(getParam<Real>("Nv")),
  _Nc(getParam<Real>("Nc")),
   _T(getParam<Real>("T")),
   _Kb(getParam<Real>("Kb")),
   _q(getParam<Real>("q")),
   _b(getParam<Real>("b"))
{
}

Real
CarrierRec::computeQpResidual()
{
  Real Rrec = 0.0;
  Rrec += -1 * _b * _test[_i][_qp] *

  (_Nc * std::exp((( _q * _u[_qp]) - _Ec) / (_Kb * _T))) * (_Nv * std::exp(( _Ev - (_q * _u[_qp])) / (_Kb * _T)));

  return Rrec;
}

Real
CarrierRec::computeQpJacobian()
{
   return -1 * _b * _test[_i][_qp] *

   ((((_Nc * _q * _phi[_j][_qp] / (_Kb * _T)) * std::exp((( _q * _u[_qp]) - _Ec) / (_Kb * _T)))*(_Nv * std::exp(( _Ev - (_q * _u[_qp])) / (_Kb * _T)))) +

   (((-1 * _Nv * _q * _phi[_j][_qp] / (_Kb * _T)) * std::exp(( _Ev - (_q * _u[_qp])) / (_Kb * _T)))*(_Nc * std::exp((( _q * _u[_qp]) - _Ec) / (_Kb * _T)))));
}
