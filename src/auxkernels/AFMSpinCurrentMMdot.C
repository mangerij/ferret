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

#include "AFMSpinCurrentMMdot.h"
registerMooseObject("FerretApp", AFMSpinCurrentMMdot);

InputParameters AFMSpinCurrentMMdot::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the AFM spin current component corresponding to the cross product of M with dM/dt");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetization");
  params.addRequiredCoupledVar("mag_z", "The z component of the magnetization");
  params.addRequiredCoupledVar("dmag_dt_x", "The time derivative of the x component of the magnetization");
  params.addRequiredCoupledVar("dmag_dt_y", "The time derivative of the y component of the magnetization");
  params.addRequiredCoupledVar("dmag_dt_z", "The time derivative of the z component of the magnetization");
  params.addRequiredParam<unsigned int>("component", "component to calculate");
  params.addParam<Real>("factor", 1.0, "scale currents");
  return params;
}


AFMSpinCurrentMMdot::AFMSpinCurrentMMdot(const InputParameters & parameters) :
  AuxKernel(parameters),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _dmag_dt_x(coupledValue("dmag_dt_x")),
   _dmag_dt_y(coupledValue("dmag_dt_y")),
   _dmag_dt_z(coupledValue("dmag_dt_z")),
   _component(getParam<unsigned int>("component")),
   _factor(getParam<Real>("factor"))
{
}

Real
AFMSpinCurrentMMdot::computeValue()
{
  if (_component == 0)
  {
   return _factor*(_dmag_dt_z[_qp]*_mag_y[_qp]-_dmag_dt_y[_qp]*_mag_z[_qp]);
  }
  else if (_component == 1)
  {
   return _factor*(-_dmag_dt_z[_qp]*_mag_x[_qp]+_dmag_dt_x[_qp]*_mag_z[_qp]);
  }
  else if (_component == 2)
  {
   return _factor*(_dmag_dt_y[_qp]*_mag_x[_qp] - _dmag_dt_x[_qp]*_mag_y[_qp]);
  }
  else
    return 0.0;
}
