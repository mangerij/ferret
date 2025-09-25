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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "HoleCurrent.h"

registerMooseObject("FerretApp", HoleCurrent);

InputParameters
HoleCurrent::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  params.addRequiredParam<Real>("Ev", "Property name of the Conduction band energy(J)");
  params.addRequiredParam<Real>("Nv", "Effective DOS of the conduction band(T=298)");
  params.addRequiredParam<Real>("T", "temperature (K)");
  params.addRequiredParam<Real>("Kb", "Boltzmann Constant (aJ/K)");
  params.addRequiredParam<Real>("q", "eV (aJ)");
  params.addRequiredParam<Real>("mup","hole mobility");
  return params;
}

HoleCurrent::HoleCurrent(const InputParameters & parameters)
  :Kernel(parameters),
  _Ev(getParam<Real>("Ev")),
  _Nv(getParam<Real>("Nv")),
  _T(getParam<Real>("T")),
  _Kb(getParam<Real>("Kb")),
  _q(getParam<Real>("q")),
  _mup(getParam<Real>("mup"))
{
}

Real
HoleCurrent::computeQpResidual()
{
  Real Rhcur = 0.0;
  Rhcur += _mup * _grad_test[_i][_qp] * (_Nv * std::exp(( _Ev - (_q * _u[_qp])) / (_Kb * _T))) * _grad_u[_qp];

  return Rhcur;
}

Real
HoleCurrent::computeQpJacobian()
{
   return  _mup * _grad_test[_i][_qp] *

   (((-1 * _Nv * _q * _phi[_j][_qp] / (_Kb * _T)) * std::exp((( _Ev - (_q * _u[_qp])) / (_Kb * _T))) * _grad_u[_qp]) +

   ((_Nv * std::exp(( _Ev - (_q * _u[_qp])) / (_Kb * _T))) * _grad_phi[_j][_qp]));
}
