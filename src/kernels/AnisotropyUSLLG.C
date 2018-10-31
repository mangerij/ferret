/**
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) a_ny later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT A_ny WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for _More details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. _Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "AnisotropyUSLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", AnisotropyUSLLG);

template<>
InputParameters validParams<AnisotropyUSLLG>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for theta, 1 for phi)");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  params.addRequiredParam<Real>("alpha", "the damping coefficient in the LLG equation");
  params.addRequiredParam<Real>("K1", "K1");
  params.addRequiredParam<Real>("K2", "K2");
  params.addRequiredParam<Real>("g0", "g0");
  params.addRequiredParam<Real>("Ms", "Ms");
  params.addRequiredParam<Real>("mu0", "mu0");
  return params;
}

AnisotropyUSLLG::AnisotropyUSLLG(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _azimuth_phi_var(coupled("azimuth_phi")),
  _polar_theta_var(coupled("polar_theta")),
  _azimuth_phi(coupledValue("azimuth_phi")),
  _polar_theta(coupledValue("polar_theta")),
  _azimuth_phi_grad(coupledGradient("azimuth_phi")),
  _polar_theta_grad(coupledGradient("polar_theta")),
  _alpha(getParam<Real>("alpha")),
  _K1(getParam<Real>("K1")),
  _K2(getParam<Real>("K2")),
  _g0(getParam<Real>("g0")),
  _Ms(getParam<Real>("Ms")),
  _mu0(getParam<Real>("mu0"))
{
}

Real
AnisotropyUSLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return (-(_g0*Utility::pow<2>(_Ms)*_test[_i][_qp]*std::sin(_polar_theta[_qp])*(8*(1 + Utility::pow<2>(_alpha))*(2*_K1 + _K2*Utility::pow<2>(_Ms) + _K2*Utility::pow<2>(_Ms)*std::cos(2*_polar_theta[_qp]))*std::sin(4*_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 
        _alpha*std::cos(_polar_theta[_qp])*(8*_K1 - _K2*Utility::pow<2>(_Ms) + 56*_K1*std::cos(2*_polar_theta[_qp]) + 4*_K2*Utility::pow<2>(_Ms)*std::cos(2*_polar_theta[_qp]) - 3*_K2*Utility::pow<2>(_Ms)*std::cos(4*_polar_theta[_qp]) - 4*std::cos(4*_azimuth_phi[_qp])*(4*_K1 + _K2*Utility::pow<2>(_Ms) + 3*_K2*Utility::pow<2>(_Ms)*std::cos(2*_polar_theta[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])))))/
   (32.*(1 + Utility::pow<2>(_alpha))*_mu0));
  }
  else if (_component == 1)
  {
    return ((_g0*Utility::pow<2>(_Ms)*_test[_i][_qp]*(-8*_alpha*(2*_K1 + _K2*Utility::pow<2>(_Ms) + _K2*Utility::pow<2>(_Ms)*std::cos(2*_polar_theta[_qp]))*std::sin(4*_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 
       (1 + Utility::pow<2>(_alpha))*std::cos(_polar_theta[_qp])*(8*_K1 - _K2*Utility::pow<2>(_Ms) + 56*_K1*std::cos(2*_polar_theta[_qp]) + 4*_K2*Utility::pow<2>(_Ms)*std::cos(2*_polar_theta[_qp]) - 3*_K2*Utility::pow<2>(_Ms)*std::cos(4*_polar_theta[_qp]) - 
          4*std::cos(4*_azimuth_phi[_qp])*(4*_K1 + _K2*Utility::pow<2>(_Ms) + 3*_K2*Utility::pow<2>(_Ms)*std::cos(2*_polar_theta[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])))))/(32.*(1 + Utility::pow<2>(_alpha))*_mu0));
  }
  else
    return 0.0;
}

Real
AnisotropyUSLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return ((_g0*Utility::pow<2>(_Ms)*_phi[_j][_qp]*_test[_i][_qp]*(-2*_alpha*(14*_K1 + _K2*Utility::pow<2>(_Ms) + (2*_K1 - _K2*Utility::pow<2>(_Ms))*std::cos(4*_azimuth_phi[_qp]))*std::cos(4*_polar_theta[_qp]) + _alpha*std::cos(2*_polar_theta[_qp])*(-8*_K1 - 5*_K2*Utility::pow<2>(_Ms) + 9*_K2*Utility::pow<2>(_Ms)*std::cos(4*_polar_theta[_qp]))*Utility::pow<2>(std::sin(2*_azimuth_phi[_qp])) - 
       4*(1 + Utility::pow<2>(_alpha))*std::cos(_polar_theta[_qp])*(6*_K1 + _K2*Utility::pow<2>(_Ms) + 5*_K2*Utility::pow<2>(_Ms)*std::cos(2*_polar_theta[_qp]))*std::sin(4*_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp]))))/(16.*(1 + Utility::pow<2>(_alpha))*_mu0));
  }
  else if (_component == 1)
  {
    return ((_g0*Utility::pow<2>(_Ms)*_phi[_j][_qp]*_test[_i][_qp]*(-2*_alpha*std::cos(4*_azimuth_phi[_qp])*(2*_K1 + _K2*Utility::pow<2>(_Ms) + _K2*Utility::pow<2>(_Ms)*std::cos(2*_polar_theta[_qp])) + (1 + Utility::pow<2>(_alpha))*std::cos(_polar_theta[_qp])*(4*_K1 + _K2*Utility::pow<2>(_Ms) + 3*_K2*Utility::pow<2>(_Ms)*std::cos(2*_polar_theta[_qp]))*std::sin(4*_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])))/
   (2.*(1 + Utility::pow<2>(_alpha))*_mu0));
  }
  else
    return 0.0;
}

Real
AnisotropyUSLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuth_phi_var)
    {
      return (-(_g0*Utility::pow<2>(_Ms)*_phi[_j][_qp]*_test[_i][_qp]*(2*(1 + Utility::pow<2>(_alpha))*std::cos(4*_azimuth_phi[_qp])*(2*_K1 + _K2*Utility::pow<2>(_Ms) + _K2*Utility::pow<2>(_Ms)*std::cos(2*_polar_theta[_qp])) + _alpha*std::cos(_polar_theta[_qp])*(4*_K1 + _K2*Utility::pow<2>(_Ms) + 3*_K2*Utility::pow<2>(_Ms)*std::cos(2*_polar_theta[_qp]))*std::sin(4*_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp])))/
   (2.*(1 + Utility::pow<2>(_alpha))*_mu0));
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
      return ((_g0*Utility::pow<2>(_Ms)*_phi[_j][_qp]*_test[_i][_qp]*((1 + Utility::pow<2>(_alpha))*(-120*_K1 + 5*_K2*Utility::pow<2>(_Ms) - (8*_K1 + 5*_K2*Utility::pow<2>(_Ms))*std::cos(4*_azimuth_phi[_qp]) - 12*(14*_K1 - _K2*Utility::pow<2>(_Ms) + (2*_K1 + _K2*Utility::pow<2>(_Ms))*std::cos(4*_azimuth_phi[_qp]))*std::cos(2*_polar_theta[_qp]) + 
          120*_K2*Utility::pow<2>(_Ms)*Utility::pow<2>(std::cos(_azimuth_phi[_qp]))*std::cos(4*_polar_theta[_qp])*Utility::pow<2>(std::sin(_azimuth_phi[_qp])))*std::sin(_polar_theta[_qp]) - 16*_alpha*(_K1 + _K2*Utility::pow<2>(_Ms)*std::cos(2*_polar_theta[_qp]))*std::sin(4*_azimuth_phi[_qp])*std::sin(2*_polar_theta[_qp])))/(32.*(1 + Utility::pow<2>(_alpha))*_mu0));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
