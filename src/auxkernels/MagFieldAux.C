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

#include "MagFieldAux.h"

registerMooseObject("FerretApp", MagFieldAux);

InputParameters MagFieldAux::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Converts polar and azimuthal solution variables to the locally saturated magnetization vector");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addCoupledVar("polar_th", "The polar angle variable");
  params.addCoupledVar("azimuthal_ph", "The azimuthal angle variable");
  return params;
}


MagFieldAux::MagFieldAux(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
  _azimuthal_ph(coupledValue("azimuthal_ph")),
  _polar_th(coupledValue("polar_th"))
{
}

Real
MagFieldAux::computeValue()
{
  if (_component == 0)
  {
   return std::sin(_polar_th[_qp]) * std::cos(_azimuthal_ph[_qp]);
  }
  else if (_component == 1)
  {
    return std::sin(_polar_th[_qp]) * std::sin(_azimuthal_ph[_qp]);
  }
  else if (_component == 2)
  {
    return std::cos(_polar_th[_qp]);
  }
  else
    return 0.0;
}
