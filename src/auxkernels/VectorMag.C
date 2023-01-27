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

#include "VectorMag.h"
#include <math.h>

registerMooseObject("FerretApp", VectorMag);

InputParameters VectorMag::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("vector_x", "The x component of the vector");
  params.addRequiredCoupledVar("vector_y", "The y component of the vector");
  params.addCoupledVar("vector_z", 0.0, "The z component of the vector");
  return params;
}


VectorMag::VectorMag(const InputParameters & parameters) :
  AuxKernel(parameters),
  _vector_x(coupledValue("vector_x")),
  _vector_y(coupledValue("vector_y")),
  _vector_z(coupledValue("vector_z"))
{}

Real
VectorMag::computeValue()
{
  RealVectorValue w(_vector_x[_qp], _vector_y[_qp], _vector_z[_qp]);
  return sqrt(w*w);
}
