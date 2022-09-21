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

#include "AFMSpinCurrentLMdot.h"
registerMooseObject("FerretApp", AFMSpinCurrentLMdot);

InputParameters AFMSpinCurrentLMdot::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the cross product of L with dm/dt");
  params.addRequiredCoupledVar("Neel_L_x", "The x component of the Neel order");
  params.addRequiredCoupledVar("Neel_L_y", "The y component of the Neel order");
  params.addRequiredCoupledVar("Neel_L_z", "The z component of the Neel order");
  params.addRequiredCoupledVar("dSSmag_dt_x", "The time derivative of the x component of the ss mag order");
  params.addRequiredCoupledVar("dSSmag_dt_y", "The time derivative of the y component of the ss mag order");
  params.addRequiredCoupledVar("dSSmag_dt_z", "The time derivative of the z component of the ss mag order");
  params.addRequiredParam<unsigned int>("component", "component to calculate");
  params.addParam<Real>("factor", 1.0, "scale currents");
  return params;
}


AFMSpinCurrentLMdot::AFMSpinCurrentLMdot(const InputParameters & parameters) :
  AuxKernel(parameters),
   _Neel_L_x(coupledValue("Neel_L_x")),
   _Neel_L_y(coupledValue("Neel_L_y")),
   _Neel_L_z(coupledValue("Neel_L_z")),
   _dSSmag_dt_x(coupledValue("dSSmag_dt_x")),
   _dSSmag_dt_y(coupledValue("dSSmag_dt_y")),
   _dSSmag_dt_z(coupledValue("dSSmag_dt_z")),
   _component(getParam<unsigned int>("component")),
   _factor(getParam<Real>("factor"))
{
}

Real
AFMSpinCurrentLMdot::computeValue()
{
  if (_component == 0)
  {
   return _factor*(_dSSmag_dt_z[_qp]*_Neel_L_y[_qp] - _dSSmag_dt_y[_qp]*_Neel_L_z[_qp]);
  }
  else if (_component == 1)
  {
   return _factor*(-_dSSmag_dt_z[_qp]*_Neel_L_x[_qp] + _dSSmag_dt_x[_qp]*_Neel_L_z[_qp]);
  }
  else if (_component == 2)
  {
   return _factor*(_dSSmag_dt_y[_qp]*_Neel_L_x[_qp] - _dSSmag_dt_x[_qp]*_Neel_L_y[_qp]);
  }
  else
    return 0.0;
}

