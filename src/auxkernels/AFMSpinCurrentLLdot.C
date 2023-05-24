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

#include "AFMSpinCurrentLLdot.h"
registerMooseObject("FerretApp", AFMSpinCurrentLLdot);

InputParameters AFMSpinCurrentLLdot::validParams()

{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the AFM spin current component corresponding to the cross product of $\\mathbf{L}$ with d$\\mathbf{L}$/dt.");
  params.addRequiredCoupledVar("Neel_L_x", "The x component of the Neel order");
  params.addRequiredCoupledVar("Neel_L_y", "The y component of the Neel order");
  params.addRequiredCoupledVar("Neel_L_z", "The z component of the Neel order");
  params.addRequiredCoupledVar("dL_dt_x", "The time derivative of the x component of the Neel order");
  params.addRequiredCoupledVar("dL_dt_y", "The time derivative of the y component of the Neel order");
  params.addRequiredCoupledVar("dL_dt_z", "The time derivative of the z component of the Neel order");
  params.addRequiredParam<unsigned int>("component", "component to calculate");
  params.addParam<Real>("factor", 1.0, "scale currents");
  return params;
}


AFMSpinCurrentLLdot::AFMSpinCurrentLLdot(const InputParameters & parameters) :
  AuxKernel(parameters),
   _Neel_L_x(coupledValue("Neel_L_x")),
   _Neel_L_y(coupledValue("Neel_L_y")),
   _Neel_L_z(coupledValue("Neel_L_z")),
   _dL_dt_x(coupledValue("dL_dt_x")),
   _dL_dt_y(coupledValue("dL_dt_y")),
   _dL_dt_z(coupledValue("dL_dt_z")),
   _component(getParam<unsigned int>("component")),
   _factor(getParam<Real>("factor"))
{
}

Real
AFMSpinCurrentLLdot::computeValue()
{
  if (_component == 0)
  {
   return _factor*(_dL_dt_z[_qp]*_Neel_L_y[_qp]-_dL_dt_y[_qp]*_Neel_L_z[_qp]);
  }
  else if (_component == 1)
  {
   return _factor*(-_dL_dt_z[_qp]*_Neel_L_x[_qp]+_dL_dt_x[_qp]*_Neel_L_z[_qp]);
  }
  else if (_component == 2)
  {
   return _factor*(_dL_dt_y[_qp]*_Neel_L_x[_qp] - _dL_dt_x[_qp]*_Neel_L_y[_qp]);
  }
  else
    return 0.0;
}
