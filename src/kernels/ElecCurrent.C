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

#include "ElecCurrent.h"

registerMooseObject("FerretApp", ElecCurrent);

InputParameters
ElecCurrent::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  params.addRequiredParam<Real>("mun", "electron mobility");
  params.addRequiredParam<Real>("Nc","Effective DOS of the conduction band(T=298)");
  params.addRequiredParam<Real>("q", "eV (aJ)");
  params.addRequiredParam<Real>("Ec", "Property name of the Conduction band energy(J)");
  params.addRequiredParam<Real>("Kb", "Boltzmann Constant (aJ/K)");
  params.addRequiredParam<Real>("T", "temperature (K)");
  return params;
}

ElecCurrent::ElecCurrent(const InputParameters & parameters)
  :Kernel(parameters),
   _mun(getParam<Real>("mun")),
   _Nc(getParam<Real>("Nc")),
   _q(getParam<Real>("q")),
   _Ec(getParam<Real>("Ec")),
   _Kb(getParam<Real>("Kb")),
   _T(getParam<Real>("T"))
{
}

Real
ElecCurrent::computeQpResidual()
{
  Real Recur = 0.0;
  Recur += -_mun * _grad_test[_i][_qp] * (_Nc * std::exp((( _q * _u[_qp]) - _Ec) / (_Kb * _T))) * _grad_u[_qp];

  return Recur;
}

Real
ElecCurrent::computeQpJacobian()
{
   return (-_mun * _grad_test[_i][_qp]) *

   (((_Nc * _q * _phi[_j][_qp] / (_Kb * _T)) * std::exp((( _q * _u[_qp]) - _Ec) / (_Kb * _T))* _grad_u[_qp]) +

   ((_Nc * std::exp(((_q * _u[_qp]) - _Ec) / (_Kb * _T))) * _grad_phi[_j][_qp]));
}
