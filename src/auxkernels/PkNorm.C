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

#include "PkNorm.h"
#include <math.h>

registerMooseObject("FerretApp", PkNorm);

template<>
InputParameters validParams<PkNorm>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("component", "the component of the normalized vector to store");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}


PkNorm::PkNorm(const InputParameters & parameters) :
  AuxKernel(parameters),
  _component(getParam<unsigned int>("component")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z"))
{}

Real
PkNorm::computeValue()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  if (_component == 0)
    return _polar_x[_qp] / sqrt(w*w);
  else if (_component == 1)
    return _polar_y[_qp] / sqrt(w*w);
  else if (_component == 2)
    return _polar_z[_qp] / sqrt(w*w);
  else
    return 0.0;
}
