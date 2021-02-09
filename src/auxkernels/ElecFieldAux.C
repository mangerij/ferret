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

#include "ElecFieldAux.h"

registerMooseObject("FerretApp", ElecFieldAux);

template<>

InputParameters validParams<ElecFieldAux>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Converts electrostatic potential to the vector electric field.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addCoupledVar("potential_E_int", "The internal electric potential variable");
  params.addCoupledVar("potential_E_ext", "The external electric potential variable");
  return params;
}


ElecFieldAux::ElecFieldAux(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _potential_E_int_grad(coupledGradient("potential_E_int")),
   _potential_E_ext_grad(coupledGradient("potential_E_ext"))
{
}

Real
ElecFieldAux::computeValue()
{
    return - _potential_E_int_grad[_qp](_component) - _potential_E_ext_grad[_qp](_component);
}
