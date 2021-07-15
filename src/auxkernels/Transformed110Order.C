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

#include "Transformed110Order.h"
#include <math.h>

registerMooseObject("FerretApp", Transformed110Order);

template<>
InputParameters validParams<Transformed110Order>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<bool>("inverse", "If this is true then inverse transformation is calculated");
  params.addRequiredParam<unsigned int>("component", "the component of the transformed vector to store");
  params.addRequiredCoupledVar("order_param_x", "The x component of the order parameter");
  params.addRequiredCoupledVar("order_param_y", "The y component of the order parameter");
  params.addCoupledVar("order_param_z", 0.0, "The z component of the order parameter");
  return params;
}


Transformed110Order::Transformed110Order(const InputParameters & parameters) :
  AuxKernel(parameters),
  _inverse(parameters.get<bool>("inverse")),
  _component(getParam<unsigned int>("component")),
  _order_param_x(coupledValue("order_param_x")),
  _order_param_y(coupledValue("order_param_y")),
  _order_param_z(coupledValue("order_param_z"))
{}

Real
Transformed110Order::computeValue()
{
//
// TODO: Note that there is no reason this needs to be hardcoded, but this will be the first step. 
//       in general, this procedure should work for any transformation
//
  if (_inverse == false)
  {
  // This takes the op_x, op_y, and op_z components back to the x,y,z orientation [110]->[100].
    if (_component == 0)
      return 0.0;
    else if (_component == 1)
      return 0.0;
    else if (_component == 2)
      return 0.0;
    else
      return 0.0;
  }
  else if (_inverse == true)
  {
  // This takes the x, y, and z components back to the op_x,op_y,op_z [110] orientation [100]->[110].
    if (_component == 0)
      return 0.0;
    else if (_component == 1)
      return 0.0;
    else if (_component == 2)
      return 0.0;
    else
      return 0.0;
  }
  else
    return 0.0;
}
