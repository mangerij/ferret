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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "HomogeneousDisplacement.h"

registerMooseObject("FerretApp", HomogeneousDisplacement);

template<>
InputParameters validParams<HomogeneousDisplacement>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("component", "the component of the normalized vector to store");
  params.addRequiredCoupledVar("disp_x", "The x component of the total elastic displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the total elastic displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the total elastic displacement");
  params.addRequiredCoupledVar("u_x", "The x component of the local elastic displacement");
  params.addRequiredCoupledVar("u_y", "The y component of the local elastic displacement");
  params.addCoupledVar("u_z", 0.0, "The z component of the local elastic displacement");
  return params;
}


HomogeneousDisplacement::HomogeneousDisplacement(const InputParameters & parameters) :
  AuxKernel(parameters),
  _component(getParam<unsigned int>("component")),
  _disp_x(coupledValue("disp_x")),
  _disp_y(coupledValue("disp_y")),
  _disp_z(coupledValue("disp_z")),
  _u_x(coupledValue("u_x")),
  _u_y(coupledValue("u_y")),
  _u_z(coupledValue("u_z"))
{
}

Real
HomogeneousDisplacement::computeValue()
{
  if (_component == 0)
    return _disp_x[_qp] - _u_x[_qp];
  else if (_component == 1)
    return _disp_y[_qp] - _u_y[_qp];
  else if (_component == 2)
    return _disp_z[_qp] - _u_z[_qp];
  else
    return 0.0;
}
