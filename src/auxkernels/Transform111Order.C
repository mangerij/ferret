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

#include "Transform111Order.h"
#include <math.h>

registerMooseObject("FerretApp", Transform111Order);

InputParameters Transform111Order::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<unsigned int>("component", "the component of the transformed vector to store");
  params.addRequiredCoupledVar("order_param_x", "The x component of the order parameter");
  params.addRequiredCoupledVar("order_param_y", "The y component of the order parameter");
  params.addCoupledVar("order_param_z", 0.0, "The z component of the order parameter");
  return params;
}


Transform111Order::Transform111Order(const InputParameters & parameters) :
  AuxKernel(parameters),
  _component(getParam<unsigned int>("component")),
  _order_param_x(coupledValue("order_param_x")),
  _order_param_y(coupledValue("order_param_y")),
  _order_param_z(coupledValue("order_param_z"))
{}

Real
Transform111Order::computeValue()
{
  if (_component == 0)
    return 0.40824829046386301637*_order_param_x[_qp] + 0.40824829046386301637*_order_param_y[_qp] - 
 0.81649658092772603273*_order_param_z[_qp];
  else if (_component == 1)
    return -0.70710678118654752440*_order_param_x[_qp] + 0.70710678118654752440*_order_param_y[_qp];
  else if (_component == 2)
    return 0.57735026918962576451*_order_param_x[_qp] + 0.57735026918962576451*_order_param_y[_qp] + 
 0.57735026918962576451*_order_param_z[_qp];
  else
    return 0.0;
}
