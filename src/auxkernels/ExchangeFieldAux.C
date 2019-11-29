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
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("magnetic_x", "The x component of the magnetization");
  params.addRequiredCoupledVar("magnetic_y", "The y component of the magnetization");
  params.addCoupledVar("magnetic_z", 0.0, "The z component of the magnetization");
  params.addRequiredParam<Real>("Ae", "Ae");
  params.addRequiredParam<Real>("Ms", "Ms");
  return params;
}


ExchangeFieldAux::ExchangeFieldAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _component(getParam<unsigned int>("component")),
  _magnetic_x_lap(coupledSecond("magnetic_x")),
  _magnetic_y_lap(coupledSecond("magnetic_y")),
  _magnetic_z_lap(coupledSecond("magnetic_z")),
  _Ae(getParam<Real>("Ae")),
  _Ms(getParam<Real>("Ms"))
{
}

Real
ExchangeFieldAux::computeValue()
{
  if (_component == 0)
    return (2.0*_Ae/_Ms)*(_magnetic_x_lap[_qp](0,0)+_magnetic_x_lap[_qp](1,1)+_magnetic_x_lap[_qp](2,2));
  else if (_component == 1)
    return (2.0*_Ae/_Ms)*(_magnetic_y_lap[_qp](0,0)+_magnetic_y_lap[_qp](1,1)+_magnetic_y_lap[_qp](2,2));
  else if (_component == 2)
    return (2.0*_Ae/_Ms)*(_magnetic_z_lap[_qp](0,0)+_magnetic_z_lap[_qp](1,1)+_magnetic_z_lap[_qp](2,2));
  else
    return 0.0;
}
