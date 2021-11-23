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

#include "InPlaneP.h"
registerMooseObject("FerretApp", InPlaneP);

InputParameters InPlaneP::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Computes some in-plane (assuming th x-y plane is the plane) polarization components");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  return params;
}

InPlaneP::InPlaneP(const InputParameters & parameters) :
  AuxKernel(parameters),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y"))
{
}

Real
InPlaneP::computeValue()
{
    return (_polar_x[_qp]*_polar_x[_qp] + _polar_y[_qp]*_polar_y[_qp]);
}
