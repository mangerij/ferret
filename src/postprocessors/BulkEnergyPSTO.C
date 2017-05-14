/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#include "BulkEnergyPSTO.h"

template<>
InputParameters validParams<BulkEnergyPSTO>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1", "The coefficients of the Landau expansion"); 
  params.addRequiredParam<Real>("alpha2", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha3", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha4", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha5", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("T", "Temperature");
  params.addRequiredParam<Real>("Tc", "Critical Temperature");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

BulkEnergyPSTO::BulkEnergyPSTO(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _alpha1(getParam<Real>("alpha1")),
  _alpha2(getParam<Real>("alpha2")),
  _alpha3(getParam<Real>("alpha3")),
  _alpha4(getParam<Real>("alpha4")),
  _alpha5(getParam<Real>("alpha5")),
  _T(getParam<Real>("T")),
  _Tc(getParam<Real>("Tc"))

{
}

Real
BulkEnergyPSTO::computeQpIntegral()
{
  return 
  _alpha1 * (_T-_Tc) * (std::pow(_polar_x[_qp], 2.0) + std::pow(_polar_y[_qp], 2.0)) + _alpha2 * (std::pow(_polar_x[_qp], 4.0) +std::pow(_polar_y[_qp], 4.0)) + _alpha3 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0) + _alpha4 * (std::pow(_polar_x[_qp], 6.0) + std::pow(_polar_y[_qp], 6.0)) + _alpha5 * (std::pow(_polar_x[_qp], 4.0) * std::pow(_polar_y[_qp], 2.0) + std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 4.0));
}
