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

#include "TimeUSLL.h"
#include "MooseVariable.h"

registerMooseObject("FerretApp", TimeUSLL);

template<>
InputParameters validParams<TimeUSLL>()
{
  InputParameters params = validParams<TimeKernel>();
  params.addClassDescription("Calculates a residual contribution for the left-hand-side of the Landau-Lifshitz equation");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for polar, 1 for azimuth)");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  return params;
}

TimeUSLL::TimeUSLL(const InputParameters & parameters) :
    TimeKernel(parameters),
  _component(getParam<unsigned int>("component")),
  _azimuth_phi_var(coupled("azimuth_phi")),
  _polar_theta_var(coupled("polar_theta")),
  _azimuth_phi(coupledValue("azimuth_phi")),
  _polar_theta(coupledValue("polar_theta")),
  _azimuth_phi_dot(coupledDot("azimuth_phi")),
  _polar_theta_dot(coupledDot("polar_theta")),
  _azimuth_phi_d_dot(coupledDotDu("azimuth_phi")),
  _polar_theta_d_dot(coupledDotDu("polar_theta"))

{
}

Real
TimeUSLL::computeQpResidual()
{
  if (_component == 0)
  {
    return -_test[_i][_qp]*_polar_theta_dot[_qp]*std::sqrt(1-std::cos(_polar_theta[_qp])*std::cos(_polar_theta[_qp]));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp]*std::sin(_polar_theta[_qp])*_azimuth_phi_dot[_qp];
  }
  else
    return 0.0;
}

Real
TimeUSLL::computeQpJacobian()
{
  if (_component == 0)
  {
    return -_test[_i][_qp]*_phi[_j][_qp]*_polar_theta_d_dot[_qp]*std::sin(_polar_theta[_qp])-(_test[_i][_qp]*_phi[_j][_qp]*_polar_theta_dot[_qp]*std::sin(_polar_theta[_qp])*cos(_polar_theta[_qp])/std::sqrt(1.0-std::cos(_polar_theta[_qp])*std::cos(_polar_theta[_qp])));
  }
  else if (_component == 1)
  {
    return std::sin(_polar_theta[_qp])*_test[_i][_qp]*_phi[_j][_qp]*_azimuth_phi_d_dot[_qp];
  }
  else
    return 0.0;
}

Real
TimeUSLL::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuth_phi_var)
    {
      return 0.0;
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _polar_theta_var)
    {
      return 2.0*std::cos(_polar_theta[_qp])*std::sin(_polar_theta[_qp])*_test[_i][_qp]*_phi[_j][_qp]*_azimuth_phi_dot[_qp];
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}

