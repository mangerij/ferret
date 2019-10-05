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

#include "MagHStrong.h"

class MagHStrong;

registerMooseObject("FerretApp", MagHStrong);

template<>
InputParameters validParams<MagHStrong>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for bound magnetic charge (div M)");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  params.addRequiredParam<Real>("Ms", "Ms");
  return params;
}

MagHStrong::MagHStrong(const InputParameters & parameters)
  :Kernel(parameters),
   _azimuth_phi_var(coupled("azimuth_phi")),
   _polar_theta_var(coupled("polar_theta")),
   _azimuth_phi(coupledValue("azimuth_phi")),
   _polar_theta(coupledValue("polar_theta")),
   _Ms(getParam<Real>("Ms"))
{
}

Real
MagHStrong::computeQpResidual()
{
  return -(_Ms*(_grad_test[_i][_qp](2)*std::cos(_polar_theta[_qp])+(_grad_test[_i][_qp](0)*std::cos(_azimuth_phi[_qp])+_grad_test[_i][_qp](1)*std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp])));
}
Real
MagHStrong::computeQpJacobian()
{
  return 0.0;
}

Real
MagHStrong::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _polar_theta_var)
  {
    return -_Ms*_phi[_j][_qp]*(-(std::cos(_polar_theta[_qp])*(_grad_test[_i][_qp](0)*std::cos(_azimuth_phi[_qp]) + _grad_test[_i][_qp](1)*std::sin(_azimuth_phi[_qp]))) + _grad_test[_i][_qp](2)*std::sin(_polar_theta[_qp]));
  }
  else if (jvar == _azimuth_phi_var)
  {
    return -_Ms*_phi[_j][_qp]*(-(_grad_test[_i][_qp](1)*std::cos(_azimuth_phi[_qp])) + _grad_test[_i][_qp](0)*std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp]);
  }
  else
    return 0.0;
}
