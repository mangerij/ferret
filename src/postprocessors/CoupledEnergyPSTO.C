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

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "CoupledEnergyPSTO.h"

template<>
InputParameters validParams<CoupledEnergyPSTO>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("x1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x2", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x3", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x4", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x5", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("x6", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("epsilon", "Constant value for strain");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

CoupledEnergyPSTO::CoupledEnergyPSTO(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _x1(getParam<Real>("x1")),
  _x2(getParam<Real>("x2")),
  _x3(getParam<Real>("x3")),
  _x4(getParam<Real>("x4")),
  _x5(getParam<Real>("x5")),
  _x6(getParam<Real>("x6")),
  _epsilon(getParam<Real>("epsilon"))

{
}

Real
CoupledEnergyPSTO::computeQpIntegral()
{
  return
(_x1 * (std::pow(_polar_x[_qp], 2.0) +std::pow(_polar_y[_qp], 2.0)) + _x2 * (std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_y[_qp], 4.0)) + _x3 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 4.0)) * _epsilon + (_x4 * (std::pow(_polar_x[_qp], 2.0) + std::pow(_polar_y[_qp], 2.0)) + _x5 * (std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_y[_qp], 4.0))+ _x6 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0)) * _epsilon * _epsilon;

}
