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

#include "ValueAux.h"
registerMooseObject("FerretApp", ValueAux);

InputParameters ValueAux::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Stores a variable as an Aux field");
  params.addRequiredCoupledVar("var1", "The first variable");
  return params;
}


ValueAux::ValueAux(const InputParameters & parameters) :
  AuxKernel(parameters),
   _var1(coupledValue("var1"))
{
}

Real
ValueAux::computeValue()
{
  return _var1[_qp];
}
