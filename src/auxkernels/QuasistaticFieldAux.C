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

#include "QuasistaticFieldAux.h"

registerMooseObject("FerretApp", QuasistaticFieldAux);

InputParameters QuasistaticFieldAux::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Converts potential to the vector field.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addCoupledVar("potential_int", "The internal potential variable");
  params.addCoupledVar("potential_ext", "The external potential variable");
  params.addParam<Real>("conversion_factor", 1.0, "conversion factor if needed");
  return params;
}


QuasistaticFieldAux::QuasistaticFieldAux(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _potential_int_grad(coupledGradient("potential_int")),
   _potential_ext_grad(coupledGradient("potential_ext")),
   _conversion_factor(getParam<Real>("conversion_factor"))
{
}

Real
QuasistaticFieldAux::computeValue()
{
    return -_conversion_factor*(_potential_int_grad[_qp](_component) - _potential_ext_grad[_qp](_component));
}
