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

   You should have received a co_polar_y[_qp] of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "DMInteractionEnergy.h"

template<>
InputParameters validParams<DMInteractionEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("antiferromag_L_x", "The x component of the antiferromagnetic vector");
  params.addCoupledVar("antiferromag_L_y", 0.0, "The y component of the antiferromagnetic vector");
  params.addCoupledVar("antiferromag_L_z", 0.0, "The z component of the antiferromagnetic vector");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt vector");
  params.addCoupledVar("antiferrodis_A_y", 0.0, "The y component of the antiferrodistortive tilt vector");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodistortive tilt vector");
  params.addRequiredParam<Real>("hD", "");
  params.addRequiredParam<Real>("chiP", "");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

DMInteractionEnergy::DMInteractionEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _antiferromag_L_x(coupledValue("antiferromag_L_x")),
  _antiferromag_L_y(coupledValue("antiferromag_L_y")),
  _antiferromag_L_z(coupledValue("antiferromag_L_z")),
  _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
  _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
  _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
  _hD(getParam<Real>("hD")),
  _chiP(getParam<Real>("chiP")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
DMInteractionEnergy::computeQpIntegral()
{
  return (-(_chiP*std::pow(_hD,2)*(std::pow(-(_antiferrodis_A_y[_qp]*_antiferromag_L_x[_qp]) + _antiferrodis_A_x[_qp]*_antiferromag_L_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp]*_antiferromag_L_x[_qp] - _antiferrodis_A_x[_qp]*_antiferromag_L_z[_qp],2) + std::pow(-(_antiferrodis_A_z[_qp]*_antiferromag_L_y[_qp]) + _antiferrodis_A_y[_qp]*_antiferromag_L_z[_qp],2)))/2.0);
}
