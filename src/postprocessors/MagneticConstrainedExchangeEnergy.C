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

#include "MagneticConstrainedExchangeEnergy.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<MagneticConstrainedExchangeEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over the magnetic exchange energy.");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  params.addRequiredParam<Real>("Ae", "Ae");
  return params;
}

MagneticConstrainedExchangeEnergy::MagneticConstrainedExchangeEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _azimuth_phi(coupledValue("azimuth_phi")),
  _polar_theta(coupledValue("polar_theta")),
  _azimuth_phi_grad(coupledGradient("azimuth_phi")),
  _polar_theta_grad(coupledGradient("polar_theta")),
  _Ae(getParam<Real>("Ae"))
{
}

Real
MagneticConstrainedExchangeEnergy::computeQpIntegral()
{
  RealVectorValue r = (Utility::pow<2>(_polar_theta_grad[_qp](0)), Utility::pow<2>(_polar_theta_grad[_qp](1)), Utility::pow<2>(_polar_theta_grad[_qp](2)));
  RealVectorValue s = (Utility::pow<2>(_azimuth_phi_grad[_qp](0)), Utility::pow<2>(_azimuth_phi_grad[_qp](1)), Utility::pow<2>(_azimuth_phi_grad[_qp](2)));
  return _Ae*(r(0)+r(1)+r(2)+(s(0)+s(1)+s(2))*Utility::pow<2>(std::sin(_polar_theta[_qp])));
}
