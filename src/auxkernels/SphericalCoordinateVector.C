/*
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published _var2y[_qp]
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

#include "SphericalCoordinateVector.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", SphericalCoordinateVector);

InputParameters SphericalCoordinateVector::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the spherical coordinates from a vector");
  params.addRequiredParam<unsigned int>("component", "Integer for which angle to calculate");
  params.addRequiredCoupledVar("var1x", "The first component of the vector");
  params.addRequiredCoupledVar("var1y", "The second component of the vector");
  params.addRequiredCoupledVar("var1z", "The third component of the vector");
  return params;
}


SphericalCoordinateVector::SphericalCoordinateVector(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _var1x(coupledValue("var1x")),
   _var1y(coupledValue("var1y")),
   _var1z(coupledValue("var1z"))
{
}

Real
SphericalCoordinateVector::computeValue()
{
  if (_component == 0) //azimuthal phi
  {
   return std::atan(_var1y[_qp]/_var1x[_qp]);
  }
  else if (_component == 1) //polar theta
  {
   return std::acos(_var1z[_qp]/std::sqrt(Utility::pow<2>(_var1x[_qp])+Utility::pow<2>(_var1y[_qp])+Utility::pow<2>(_var1z[_qp])));
  }
  else
    return 0.0;
}
