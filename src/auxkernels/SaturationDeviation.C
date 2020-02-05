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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "SaturationDeviation.h"
#include <math.h>

registerMooseObject("FerretApp", SaturationDeviation);

template<>
InputParameters validParams<SaturationDeviation>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the magnetization");
  params.addRequiredParam<Real>("Ms", "the saturation magnetization");
  return params;
}


SaturationDeviation::SaturationDeviation(const InputParameters & parameters) :
  AuxKernel(parameters),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _Ms(getParam<Real>("Ms"))
{}

Real
SaturationDeviation::computeValue()
{
  RealVectorValue w(_mag_x[_qp], _mag_y[_qp], _mag_z[_qp]);
  return sqrt(w*w) - _Ms;
}
