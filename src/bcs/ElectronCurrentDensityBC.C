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

#include "ElectronCurrentDensityBC.h"

registerMooseObject("FerretApp", ElectronCurrentDensityBC);

InputParameters ElectronCurrentDensityBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredParam<Real>("Ec", "Property name of the Conduction band energy(J)");
  params.addRequiredParam<Real>("Nc", "Effective DOS of the conduction band(T=298)");
  params.addRequiredParam<Real>("T", "temperature (K)");
  params.addRequiredParam<Real>("Kb", "Boltzmann Constant (aJ/K)");
  params.addRequiredParam<Real>("q", "eV (aJ)");
  params.addRequiredParam<Real>("mun","hole mobility");
  return params;
}

ElectronCurrentDensityBC::ElectronCurrentDensityBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _Ec(getParam<Real>("Ec")),
  _Nc(getParam<Real>("Nc")),
  _T(getParam<Real>("T")),
  _Kb(getParam<Real>("Kb")),
  _q(getParam<Real>("q")),
  _mun(getParam<Real>("mun"))
{}

Real
ElectronCurrentDensityBC::computeQpResidual()
{
  Real ECDBC = 0;

  for(int i = 0; i < 3; ++i){
    ECDBC += -_q * _mun * _Nc * std::exp(((_q * _u[_qp]) - _Ec) / (_Kb * _T)) * _grad_u[_qp](i) -

    ((_Kb * _T)/_q) * _Nc * std::exp(((_q * _grad_u[_qp](i)) - _Ec) / (_Kb * _T)) * _normals[_qp](i);
  }
  return ECDBC;
}
