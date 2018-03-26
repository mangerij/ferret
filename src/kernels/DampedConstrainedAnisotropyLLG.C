/**
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for _More details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. _Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "DampedConstrainedAnisotropyLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", DampedConstrainedAnisotropyLLG);

template<>
InputParameters validParams<DampedConstrainedAnisotropyLLG>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  params.addRequiredParam<Real>("alpha", "the damping coefficient in the LLG equation");
  params.addRequiredParam<Real>("K1", "K1");
  params.addRequiredParam<Real>("K2", "K2");
  params.addRequiredParam<Real>("g0", "g0");
  params.addRequiredParam<Real>("Ae", "Ae");
  params.addRequiredParam<Real>("M", "M");
  return params;
}

DampedConstrainedAnisotropyLLG::DampedConstrainedAnisotropyLLG(const InputParameters & parameters)
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
  _Ae(getParam<Real>("Ae")),
  _M(getParam<Real>("M"))
{
}

Real
DampedConstrainedAnisotropyLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * ((-0.5*_g0*(4.*_alpha*_K1*Utility::pow<3>(std::cos(_polar_theta[_qp]))*std::sin(_polar_theta[_qp]) + 8.*_K2*std::cos(_azimuth_phi[_qp])*Utility::pow<2>(std::cos(_polar_theta[_qp]))*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp])) + _K1*std::sin(4.*_azimuth_phi[_qp])*Utility::pow<3>(std::sin(_polar_theta[_qp])) - 
       1.*_alpha*std::cos(_polar_theta[_qp])*std::sin(_polar_theta[_qp])*(-4.*_Ae*(Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2))) + _K1*(3. + std::cos(4.*_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 4.*_K2*Utility::pow<4>(std::sin(_azimuth_phi[_qp]))*Utility::pow<4>(std::sin(_polar_theta[_qp]))) + 2.*_Ae*_alpha*Utility::pow<2>(_azimuth_phi_grad[_qp](1))*std::sin(2.*_polar_theta[_qp]) + 
       _alpha*_K2*Utility::pow<4>(std::sin(_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(2.*_polar_theta[_qp]))))/((1. + Utility::pow<2>(_alpha))*_M));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * ((0.5*_g0*(_K1*std::cos(3.*_polar_theta[_qp]) - 8.*_alpha*_K2*std::cos(_azimuth_phi[_qp])*Utility::pow<2>(std::cos(_polar_theta[_qp]))*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) - 1.*_alpha*_K1*std::sin(4.*_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 
       std::cos(_polar_theta[_qp])*(4.*_Ae*(Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2))) + 3.*_K1 - 1.*_K1*(3. + std::cos(4.*_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) - 4.*_K2*Utility::pow<4>(std::sin(_azimuth_phi[_qp]))*Utility::pow<4>(std::sin(_polar_theta[_qp]))) + 
       _K2*(1.0/std::sin(_polar_theta[_qp]))*Utility::pow<4>(std::sin(_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(2.*_polar_theta[_qp]))))/((1. + Utility::pow<2>(_alpha))*_M));
  }
  else
    return 0.0;
}

Real
DampedConstrainedAnisotropyLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp]*((-0.5*_g0*_phi[_j][_qp]*(4.*_alpha*_K1*Utility::pow<4>(std::cos(_polar_theta[_qp])) + 24.*_K2*std::cos(_azimuth_phi[_qp])*Utility::pow<3>(std::cos(_polar_theta[_qp]))*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 
       std::cos(_polar_theta[_qp])*(3.*_K1*std::sin(4.*_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp])) - 16.*_K2*std::cos(_azimuth_phi[_qp])*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*Utility::pow<4>(std::sin(_polar_theta[_qp]))) + 
       _alpha*Utility::pow<2>(std::cos(_polar_theta[_qp]))*(4.*_Ae*(Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2))) - 3.*_K1*(7. + std::cos(4.*_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) - 20.*_K2*Utility::pow<4>(std::sin(_azimuth_phi[_qp]))*Utility::pow<4>(std::sin(_polar_theta[_qp]))) + 
       _alpha*(4.*_Ae*Utility::pow<2>(_azimuth_phi_grad[_qp](1))*std::cos(2.*_polar_theta[_qp]) - 4.*_Ae*(Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2)))*Utility::pow<2>(std::sin(_polar_theta[_qp])) + _K1*(3. + std::cos(4.*_azimuth_phi[_qp]))*Utility::pow<4>(std::sin(_polar_theta[_qp])) + 4.*_K2*Utility::pow<4>(std::sin(_azimuth_phi[_qp]))*Utility::pow<6>(std::sin(_polar_theta[_qp])) + 
          3.*_K2*Utility::pow<4>(std::sin(_azimuth_phi[_qp]))*std::sin(2.*_polar_theta[_qp])*std::sin(4.*_polar_theta[_qp]))))/((1. + Utility::pow<2>(_alpha))*_M));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp]*((_g0*(-2.*_alpha*_K1*_phi[_j][_qp]*std::cos(4.*_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 4.*std::cos(_polar_theta[_qp])*(_Ae*(_azimuth_phi_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_phi[_j][_qp](2)) + _K2*_phi[_j][_qp]*std::cos(_azimuth_phi[_qp])*(1. + 3.*std::cos(2.*_polar_theta[_qp]))*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp]))) + 
       _phi[_j][_qp]*std::sin(2.*_polar_theta[_qp])*(_K1*std::sin(4.*_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]) - 1.*_alpha*_K2*(1. + 2.*std::cos(2.*_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*std::sin(2.*_polar_theta[_qp]))))/((1. + Utility::pow<2>(_alpha))*_M));
  }
  else
    return 0.0;
}

Real
DampedConstrainedAnisotropyLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuth_phi_var)
    {
      return _test[_i][_qp]*((-2*_g0*_phi[_j][_qp]*(_K1*std::cos(4*_azimuth_phi[_qp]) + _K2*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*(6*Utility::pow<2>(std::cos(_azimuth_phi[_qp]))*Utility::pow<2>(std::cos(_polar_theta[_qp])) + _alpha*std::cos(_azimuth_phi[_qp])*(5*std::cos(_polar_theta[_qp]) + 3*std::cos(3*_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) - 2*Utility::pow<2>(std::cos(_polar_theta[_qp]))*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))) + _alpha*_K1*std::cos(_polar_theta[_qp])*std::sin(4*_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp])) - 2*_Ae*_alpha*(_azimuth_phi_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_phi[_j][_qp](2))*_g0*std::sin(2*_polar_theta[_qp]))/((1 + Utility::pow<2>(_alpha))*_M));
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
      return _test[_i][_qp] * ((0.5*_g0*_phi[_j][_qp]*(-2.*(2.*_Ae*(Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2))) + _K1 + _K1*(7. + std::cos(4.*_azimuth_phi[_qp]))*Utility::pow<2>(std::cos(_polar_theta[_qp])) + _K1*std::cos(2.*_polar_theta[_qp]) + 8.*_alpha*_K2*std::cos(_azimuth_phi[_qp])*Utility::pow<3>(std::cos(_polar_theta[_qp]))*Utility::pow<3>(std::sin(_azimuth_phi[_qp])) - 
          8.*_K2*Utility::pow<4>(std::cos(_polar_theta[_qp]))*Utility::pow<4>(std::sin(_azimuth_phi[_qp])))*std::sin(_polar_theta[_qp]) + (3.*_K1 + _K1*std::cos(4.*_azimuth_phi[_qp]) + 8.*_K2*std::cos(_polar_theta[_qp])*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*(2.*_alpha*std::cos(_azimuth_phi[_qp]) - 5.*std::cos(_polar_theta[_qp])*std::sin(_azimuth_phi[_qp])))*Utility::pow<3>(std::sin(_polar_theta[_qp])) + 
       4.*_K2*Utility::pow<4>(std::sin(_azimuth_phi[_qp]))*Utility::pow<5>(std::sin(_polar_theta[_qp])) - 1.*_alpha*_K1*std::sin(4.*_azimuth_phi[_qp])*std::sin(2.*_polar_theta[_qp])))/((1. + Utility::pow<2>(_alpha))*_M));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
