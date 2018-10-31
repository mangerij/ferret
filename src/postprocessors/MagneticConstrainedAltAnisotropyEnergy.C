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

#include "MagneticConstrainedAltAnisotropyEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MagneticConstrainedAltAnisotropyEnergy);

template<>
InputParameters validParams<MagneticConstrainedAltAnisotropyEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over the magnetic anisotropy energy.");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  params.addRequiredParam<Real>("K1", "K1");
  params.addRequiredParam<Real>("K2", "K2");
  params.addRequiredParam<Real>("Ms", "Ms");
  return params;
}

MagneticConstrainedAltAnisotropyEnergy::MagneticConstrainedAltAnisotropyEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _azimuth_phi(coupledValue("azimuth_phi")),
  _polar_theta(coupledValue("polar_theta")),
  _K1(getParam<Real>("K1")),
  _K2(getParam<Real>("K2")),
  _Ms(getParam<Real>("Ms"))
{
}

Real
MagneticConstrainedAltAnisotropyEnergy::computeQpIntegral()
{
  return _K2*Utility::pow<6>(_Ms)*Utility::pow<2>(std::cos(_azimuth_phi[_qp]))*Utility::pow<2>(std::cos(_polar_theta[_qp]))*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*Utility::pow<4>(std::sin(_polar_theta[_qp])) + _K1*(Utility::pow<2>(_Ms)*Utility::pow<2>(std::cos(_polar_theta[_qp])) + Utility::pow<2>(_Ms)*Utility::pow<2>(std::cos(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) + Utility::pow<2>(_Ms)*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])));
}
