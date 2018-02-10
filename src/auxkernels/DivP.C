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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "DivP.h"
template<>

InputParameters validParams<DivP>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculates div P");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}


DivP::DivP(const InputParameters & parameters) :
  AuxKernel(parameters),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _polar_x_grad(coupledGradient("polar_x")),
   _polar_y_grad(coupledGradient("polar_y")),
   _polar_z_grad(coupledGradient("polar_z"))
{
}

Real
DivP::computeValue()

{
    return _polar_x_grad[_qp](0) + _polar_y_grad[_qp](1) + _polar_z_grad[_qp](2);
}
