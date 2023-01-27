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

#include "ElecGen.h"

registerMooseObject("FerretApp", ElecGen);

template<>
InputParameters validParams<ElecGen>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  params.addRequiredParam<Real>("Ec", "Property name of the Conduction band energy(J)");
  params.addRequiredParam<Real>("Nc","Effective DOS of the conduction band(T=298)");
  params.addRequiredParam<Real>("T", "temperature (K)");
  params.addRequiredParam<Real>("Kb", "Boltzmann Constant (aJ/K)");
  params.addRequiredParam<Real>("q", "eV (aJ)");
  params.addRequiredParam<Real>("mun", "electron mobility");
  return params;
}

ElecGen::ElecGen(const InputParameters & parameters)
  :Kernel(parameters),
  _Ec(getParam<Real>("Ec")),
  _Nc(getParam<Real>("Nc")),
  _T(getParam<Real>("T")),
  _Kb(getParam<Real>("Kb")),
  _q(getParam<Real>("q")),
  _mun(getParam<Real>("mun"))
{
}

Real
ElecGen::computeQpResidual()
{
  Real Regen = 0.0;
  Regen += ((_mun * _Kb * _T / _q)) *

  ((_grad_test[_i][_qp](0) * _Nc * std::exp(( (_q * _grad_u[_qp](0)) - _Ec) / (_Kb * _T))) +

  (_grad_test[_i][_qp](1) * _Nc * std::exp(((_q * _grad_u[_qp](1)) - _Ec) / (_Kb * _T))) +

  (_grad_test[_i][_qp](2) * _Nc * std::exp(((_q * _grad_u[_qp](2)) - _Ec) / (_Kb * _T))));

  return Regen;
}

Real
ElecGen::computeQpJacobian()
{
   return (( 1 *_mun * _Kb * _T / _q)) *

   (((_grad_test[_i][_qp](0) * _Nc * _q * _grad_phi[_j][_qp](0) / (_Kb * _T)) * std::exp(((_q * _grad_u[_qp](0)) - _Ec) / (_Kb * _T)))

   + ((_grad_test[_i][_qp](1) * _Nc * _q * _grad_phi[_j][_qp](1) / (_Kb * _T)) * std::exp(((_q * _grad_u[_qp](1)) - _Ec) / (_Kb * _T)))

   + ((_grad_test[_i][_qp](2) * _Nc * _q * _grad_phi[_j][_qp](2) / (_Kb * _T)) * std::exp(((_q * _grad_u[_qp](2)) - _Ec) / (_Kb * _T))));
}
