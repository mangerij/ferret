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
  params.addRequiredParam<Real>("nx", "nx");
  params.addRequiredParam<Real>("ny", "ny");
  params.addRequiredParam<Real>("nz", "nz");
  params.addRequiredParam<Real>("Ku", "Ku");
  params.addRequiredParam<Real>("M", "M");
  return params;
}

MagneticConstrainedAltAnisotropyEnergy::MagneticConstrainedAltAnisotropyEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _azimuth_phi(coupledValue("azimuth_phi")),
  _polar_theta(coupledValue("polar_theta")),
  _nx(getParam<Real>("nx")),
  _ny(getParam<Real>("ny")),
  _nz(getParam<Real>("nz")),
  _Ku(getParam<Real>("Ku")),
  _M(getParam<Real>("M"))
{
}

Real
MagneticConstrainedAltAnisotropyEnergy::computeQpIntegral()
{
  return _Ku*(Utility::pow<2>(_M)*Utility::pow<2>(_nz*std::cos(_polar_theta[_qp]) + (_nx*std::cos(_azimuth_phi[_qp]) + _ny*std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp])));
}
