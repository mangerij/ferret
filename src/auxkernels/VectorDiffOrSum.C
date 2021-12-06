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

#include "VectorDiffOrSum.h"
registerMooseObject("FerretApp", VectorDiffOrSum);

template<>

InputParameters validParams<VectorDiffOrSum>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculates the difference or sum of a variable");
  params.addRequiredCoupledVar("var1", "The first variable");
  params.addRequiredCoupledVar("var2", "The second variable");
  params.addRequiredParam<unsigned int>("diffOrSum", "A flag for diff or sum");
  return params;
}


VectorDiffOrSum::VectorDiffOrSum(const InputParameters & parameters) :
  AuxKernel(parameters),
   _var1(coupledValue("var1")),
   _var2(coupledValue("var2")),
   _diffOrSum(getParam<unsigned int>("diffOrSum"))
{
}

Real
VectorDiffOrSum::computeValue()
{
  if (_diffOrSum == 0)
  {
   return _var1[_qp]-_var2[_qp];
  }
  else if (_diffOrSum == 1)
  {
   return _var1[_qp]+_var2[_qp];
  }
  else
    return 0.0;
}
