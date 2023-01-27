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

#include "HydrostaticBC.h"

registerMooseObject("FerretApp", HydrostaticBC);

InputParameters HydrostaticBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredParam<Real>("pressure","Specify the hydrostatic pressure");
  params.addRequiredParam<int>("component","Component of displacement for BC");

  return params;
}

HydrostaticBC::HydrostaticBC(const InputParameters & parameters) :
  IntegratedBC(parameters),
  _pressure(getParam<Real>("pressure")),
  _component(getParam<int>("component"))
{}

Real
HydrostaticBC::computeQpResidual()
{
  Real traction[3];
  for(int i = 0; i < 3; ++i){
    traction[i]=-_pressure * _normals[_qp](i);
  }
  return -_test[_i][_qp] * traction[_component];
}
