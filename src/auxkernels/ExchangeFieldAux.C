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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ExchangeFieldAux.h"

registerMooseObject("FerretApp", ExchangeFieldAux);

template<>

InputParameters validParams<ExchangeFieldAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Computes the exchange field");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the magnetization");
  return params;
}


ExchangeFieldAux::ExchangeFieldAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _mag_x_grad(coupledGradient("mag_x")),
  _mag_y_grad(coupledGradient("mag_y")),
  _mag_z_grad(coupledGradient("mag_z")),
  _Ae(getMaterialProperty<Real>("Ae"))
{
}

Real
ExchangeFieldAux::computeValue()
{
  return (_Ae[_qp])*(Utility::pow<2>(_mag_x_grad[_qp](0))+Utility::pow<2>(_mag_x_grad[_qp](1))+Utility::pow<2>(_mag_x_grad[_qp](2))+Utility::pow<2>(_mag_y_grad[_qp](0))+Utility::pow<2>(_mag_y_grad[_qp](1))+Utility::pow<2>(_mag_y_grad[_qp](2))+Utility::pow<2>(_mag_z_grad[_qp](0))+Utility::pow<2>(_mag_z_grad[_qp](1))+Utility::pow<2>(_mag_z_grad[_qp](2)));
}
