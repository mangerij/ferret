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

#include "DMInteractionEnergy.h"

registerMooseObject("FerretApp", DMInteractionEnergy);

template<>
InputParameters validParams<DMInteractionEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over the DM interaction free energy density (coupling AFD and magnetic ordering).");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetization vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetization vector");
  params.addCoupledVar("mag_z", 0.0, "The z component of the magnetization vector");
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
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
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
  return (-(_chiP*std::pow(_hD,2)*(std::pow(-(_antiferrodis_A_y[_qp]*_mag_x[_qp]) + _antiferrodis_A_x[_qp]*_mag_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp]*_mag_x[_qp] - _antiferrodis_A_x[_qp]*_mag_z[_qp],2) + std::pow(-(_antiferrodis_A_z[_qp]*_mag_y[_qp]) + _antiferrodis_A_y[_qp]*_mag_z[_qp],2)))/2.0);
}
