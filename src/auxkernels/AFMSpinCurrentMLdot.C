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

#include "AFMSpinCurrentMLdot.h"
registerMooseObject("FerretApp", AFMSpinCurrentMLdot);

InputParameters AFMSpinCurrentMLdot::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the cross product of m with dL/dt");
  params.addRequiredCoupledVar("dL_dt_x", "The time derivative of the x component of the Neel order");
  params.addRequiredCoupledVar("dL_dt_y", "The time derivative of the y component of the Neel order");
  params.addRequiredCoupledVar("dL_dt_z", "The time derivative of the z component of the Neel order");
  params.addRequiredCoupledVar("SSmag_x", "The x component of the ss mag order");
  params.addRequiredCoupledVar("SSmag_y", "The y component of the ss mag order");
  params.addRequiredCoupledVar("SSmag_z", "The z component of the ss mag order");
  params.addRequiredParam<unsigned int>("component", "component to calculate");
  params.addParam<Real>("factor", 1.0, "scale currents");
  return params;
}


AFMSpinCurrentMLdot::AFMSpinCurrentMLdot(const InputParameters & parameters) :
  AuxKernel(parameters),
   _dL_dt_x(coupledValue("dL_dt_x")),
   _dL_dt_y(coupledValue("dL_dt_y")),
   _dL_dt_z(coupledValue("dL_dt_z")),
   _SSmag_x(coupledValue("SSmag_x")),
   _SSmag_y(coupledValue("SSmag_y")),
   _SSmag_z(coupledValue("SSmag_z")),
   _component(getParam<unsigned int>("component")),
   _factor(getParam<Real>("factor"))
{
}

Real
AFMSpinCurrentMLdot::computeValue()
{
  if (_component == 0)
  {
   return _factor*(_dL_dt_z[_qp]*_SSmag_y[_qp] - _dL_dt_y[_qp]*_SSmag_z[_qp]);
  }
  else if (_component == 1)
  {
   return _factor*(-_dL_dt_z[_qp]*_SSmag_x[_qp] + _dL_dt_x[_qp]*_SSmag_z[_qp]);
  }
  else if (_component == 2)
  {
   return _factor*(_dL_dt_y[_qp]*_SSmag_x[_qp] - _dL_dt_x[_qp]*_SSmag_y[_qp]);
  }
  else
    return 0.0;
}

