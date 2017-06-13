/**
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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "PyFieldAux.h"

template<>

InputParameters validParams<PyFieldAux>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<Real>("permittivity_int", "internal permittivity");
  params.addRequiredParam<Real>("permittivity_ext", "external permittivity");
  ///params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
  params.addRequiredCoupledVar("potential_ext", "The external electric potential variable");
  return params;
}


PyFieldAux::PyFieldAux(const InputParameters & parameters) :
  AuxKernel(parameters ),
   _permittivity_int(getParam<Real>("permittivity_int")),
   _permittivity_ext(getParam<Real>("permittivity_ext")),
 /// _potential_int_grad(coupledGradient("potential_int")),
   _potential_ext_grad(coupledGradient("potential_ext"))
{
}

Real
PyFieldAux::computeValue()

{
    return -( _permittivity_int - _permittivity_ext) * _potential_ext_grad[_qp](1);
}
