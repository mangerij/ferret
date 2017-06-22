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

#include "RenormalizedBulkEnergy.h"

template<>
InputParameters validParams<RenormalizedBulkEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("epsilon",0.0,"the strain coupling");
  params.addParam<Real>("T",0.0,"the temperature");
  return params;
}

RenormalizedBulkEnergy::RenormalizedBulkEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _epsilon(getParam<Real>("epsilon")),
  _T(getParam<Real>("T"))
{
}

Real
RenormalizedBulkEnergy::computeQpIntegral()
{
  return - ( 0.488524 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0) + 0.77016 * (std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_y[_qp], 4.0)) - 3.7 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) + 0.132318 * std::pow(_polar_z[_qp], 4.0) + 0.253636 * (std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) + std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0)) + 0.26 * ( std::pow(_polar_x[_qp], 6.0) + std::pow(_polar_y[_qp], 6.0) + std::pow(_polar_z[_qp], 6.0)) + 0.61 * ((std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_y[_qp], 4.0)) *std::pow(_polar_z[_qp], 2.0) + (std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_z[_qp], 4.0)) *std::pow(_polar_y[_qp], 2.0) + (std::pow(_polar_y[_qp], 4.0) + std::pow(_polar_z[_qp], 4.0)) *std::pow(_polar_x[_qp], 2.0)) + (std::pow(_polar_x[_qp], 2.0) + std::pow(_polar_y[_qp], 2.0)) * (-4.92223 * (-765.1 + _T) - 19.0909 * _epsilon) + std::pow(_polar_z[_qp], 2.0) * (-4.92223 * (-765.1 + _T) + 15.7576 * _epsilon));
  //0.488524 Px^2 Py^2 + 0.77016 (Px^4 + Py^4) - 3.7 Px^2 Py^2 Pz^2 +  0.132318 Pz^4 + 0.253636 (Px^2 Pz^2 + Py^2 Pz^2) +  0.26 (Px^6 + Py^6 + Pz^6) +  0.61 ((Px^4 + Py^4) Pz^2 + Py^2 (Px^4 + Pz^4) +     Px^2 (Py^4 + Pz^4)) + (Px^2 + Py^2) (0.000229917 -     19.0909 \[Epsilon]) + 303.03 \[Epsilon]^2 +  Pz^2 (0.000229917 + 15.7576 \[Epsilon])
}
