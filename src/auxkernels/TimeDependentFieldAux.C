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

#include "TimeDependentFieldAux.h"

registerMooseObject("FerretApp", TimeDependentFieldAux);

InputParameters TimeDependentFieldAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Adds time-dependence to a spatial-varying field");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the component of vector field to compute.");
  params.addRequiredCoupledVar("H0ext_x", "The fixed x-component of the vector field variable");
  params.addRequiredCoupledVar("H0ext_y", "The fixed y-component of the vector field variable");
  params.addRequiredCoupledVar("H0ext_z", "The fixed z-component of the vector field variable");
  params.addRequiredCoupledVar("time_dependence", "The time dependence of this field");
  return params;
}


TimeDependentFieldAux::TimeDependentFieldAux(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _H0ext_x(coupledValue("H0ext_x")),
   _H0ext_y(coupledValue("H0ext_y")),
   _H0ext_z(coupledValue("H0ext_z")),
   _time_dependence(coupledValue("time_dependence"))
{
}

Real
TimeDependentFieldAux::computeValue()
{
  if (_component == 0)
  {
    return _H0ext_x[_qp]*_time_dependence[_qp];
  }
  else if (_component == 1)
  {
    return _H0ext_y[_qp]*_time_dependence[_qp];
  }
  else if (_component == 2)
  {
    return _H0ext_z[_qp]*_time_dependence[_qp];
  }
  else
   return 0.0;
}


