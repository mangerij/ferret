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

#include "DampedConstrainedAltAnisotropyLLG.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<DampedConstrainedAltAnisotropyLLG>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  params.addRequiredParam<Real>("alpha", "the damping coefficient in the LLG equation");
  params.addRequiredParam<Real>("nx", "nx");
  params.addRequiredParam<Real>("ny", "ny");
  params.addRequiredParam<Real>("nz", "nz");
  params.addRequiredParam<Real>("Ku", "Ku");
  params.addRequiredParam<Real>("g0", "g0");
  params.addRequiredParam<Real>("Ae", "Ae");
  params.addRequiredParam<Real>("M", "M");
  return params;
}

DampedConstrainedAltAnisotropyLLG::DampedConstrainedAltAnisotropyLLG(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _azimuth_phi_var(coupled("azimuth_phi")),
  _polar_theta_var(coupled("polar_theta")),
  _azimuth_phi(coupledValue("azimuth_phi")),
  _polar_theta(coupledValue("polar_theta")),
  _azimuth_phi_grad(coupledGradient("azimuth_phi")),
  _polar_theta_grad(coupledGradient("polar_theta")),
  _alpha(getParam<Real>("alpha")),
  _nx(getParam<Real>("nx")),
  _ny(getParam<Real>("ny")),
  _nz(getParam<Real>("nz")),
  _Ku(getParam<Real>("Ku")),
  _g0(getParam<Real>("g0")),
  _Ae(getParam<Real>("Ae")),
  _M(getParam<Real>("M"))
{
}

Real
DampedConstrainedAltAnisotropyLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * ((_g0*(-2.*_Ku*Utility::pow<2>(_M)*(1.0/std::sin(_polar_theta[_qp]))*(_nz*std::cos(_azimuth_phi[_qp]) + std::sin(_azimuth_phi[_qp])*(_nx*std::cos(_polar_theta[_qp]) + _ny*std::sin(_polar_theta[_qp])))*(std::cos(_azimuth_phi[_qp])*(_nx*std::cos(_polar_theta[_qp]) + _ny*std::sin(_polar_theta[_qp])) - 1.*std::sin(_azimuth_phi[_qp])*(_nz + _alpha*std::sin(_polar_theta[_qp])*(-1.*_ny*std::cos(_polar_theta[_qp]) + _nx*std::sin(_polar_theta[_qp])))) - 1.*_Ae*_alpha*(Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2)))*std::sin(2.*_polar_theta[_qp])))/((1. + Utility::pow<2>(_alpha))*_M));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * ((_g0*(2.*std::cos(_polar_theta[_qp])*(_Ae*(Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2))) + _Ku*Utility::pow<2>(_M)*(-1.*Utility::pow<2>(_nx) + Utility::pow<2>(_ny))*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))) + 
       _Ku*Utility::pow<2>(_M)*(-2.*_alpha*_nz*std::cos(2.*_azimuth_phi[_qp])*(_ny + _nx*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])))*(1.0/std::sin(_polar_theta[_qp])) + 2.*_nx*_ny*std::sin(_azimuth_phi[_qp])*(-2.*_alpha*std::cos(_azimuth_phi[_qp])*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])) + std::cos(2.*_polar_theta[_qp])*(1.0/std::sin(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])) - 
          1.*(_alpha*Utility::pow<2>(_ny) + _nx*_nz - 1.*_ny*_nz*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])) + _alpha*Utility::pow<2>(_nx)*Utility::pow<2>((std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))) - 1.*_alpha*Utility::pow<2>(_nz)*Utility::pow<2>((1.0/std::sin(_polar_theta[_qp]))))*std::sin(2.*_azimuth_phi[_qp]))))/((1. + Utility::pow<2>(_alpha))*_M));
  }
  else
    return 0.0;
}

Real
DampedConstrainedAltAnisotropyLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * ((2.*_g0*_phi[_j][_qp]*(-1.*_Ae*_alpha*(Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2)))*std::cos(2.*_polar_theta[_qp]) + 
       _Ku*Utility::pow<2>(_M)*(_nx*_nz*Utility::pow<2>(std::cos(_azimuth_phi[_qp]))*Utility::pow<2>((1.0/std::sin(_polar_theta[_qp]))) + std::cos(_azimuth_phi[_qp])*std::sin(_azimuth_phi[_qp])*((Utility::pow<2>(_nx) - 1.*Utility::pow<2>(_ny) + _alpha*_nx*_nz)*std::cos(_polar_theta[_qp]) + (_nx - 1.*_nz)*(_nx + _nz)*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*(1.0/std::sin(_polar_theta[_qp])) + _alpha*_ny*_nz*std::sin(_polar_theta[_qp])) + std::sin(_polar_theta[_qp])*(_nx*_ny*std::sin(2.*_azimuth_phi[_qp]) + (_alpha*(_nx + _ny + (_nx - 1.*_ny)*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])))*(-1.*_nx + _ny + (_nx + _ny)*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))) - 1.*_nx*_nz*Utility::pow<4>((1.0/std::sin(_polar_theta[_qp]))))*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp])))))/((1. + Utility::pow<2>(_alpha))*_M));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * ((2.*_g0*(2.*_Ae*(_azimuth_phi_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_phi[_j][_qp](2))*std::cos(_polar_theta[_qp]) + _Ku*Utility::pow<2>(_M)*_phi[_j][_qp]*
        (-1.*std::cos(2.*_azimuth_phi[_qp])*(_alpha*Utility::pow<2>(_ny) + _nx*_nz + (std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*(2.*_alpha*_nx*_ny - 1.*_ny*_nz + _alpha*Utility::pow<2>(_nx)*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))) - 1.*_alpha*Utility::pow<2>(_nz)*Utility::pow<2>((1.0/std::sin(_polar_theta[_qp])))) + 
          (_ny + _nx*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])))*(-1.*_nx + _ny*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])) + 2.*_alpha*_nz*Utility::pow<2>((1.0/std::sin(_polar_theta[_qp]))))*std::sin(2.*_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]))))/((1. + Utility::pow<2>(_alpha))*_M));
  }
  else
    return 0.0;
}

Real
DampedConstrainedAltAnisotropyLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuth_phi_var)
    {
      return _test[_i][_qp] * ((2.*_g0*(-1.*_Ku*Utility::pow<2>(_M)*_phi[_j][_qp]*Utility::pow<2>(std::cos(_azimuth_phi[_qp]))*(Utility::pow<2>(_ny) - 1.*_alpha*_nx*_nz + (std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*(2.*_nx*_ny + _alpha*_ny*_nz + Utility::pow<2>(_nx)*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))) - 1.*Utility::pow<2>(_nz)*Utility::pow<2>((1.0/std::sin(_polar_theta[_qp]))))*std::sin(_polar_theta[_qp]) + 
       _Ku*Utility::pow<2>(_M)*_phi[_j][_qp]*((1.0/std::sin(_polar_theta[_qp]))*(Utility::pow<2>(_ny) - 1.*_alpha*_nx*_nz + (std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*(2.*_nx*_ny + _alpha*_ny*_nz + Utility::pow<2>(_nx)*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))) - 1.*Utility::pow<2>(_nz)*Utility::pow<2>((1.0/std::sin(_polar_theta[_qp]))))*Utility::pow<2>(std::sin(_azimuth_phi[_qp])) + _alpha*_nx*_ny*std::sin(2.*_azimuth_phi[_qp]))*
        Utility::pow<2>(std::sin(_polar_theta[_qp])) + 2.*_Ku*Utility::pow<2>(_M)*_phi[_j][_qp]*std::cos(_azimuth_phi[_qp])*std::sin(_azimuth_phi[_qp])*(-1.*_alpha*_nx*_ny*Utility::pow<2>(std::cos(_polar_theta[_qp])) + 2.*_nz*(_ny + _nx*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))) + _alpha*(_nx - 1.*_ny)*(_nx + _ny)*std::cos(_polar_theta[_qp])*std::sin(_polar_theta[_qp])) - 
       1.*_Ae*_alpha*(_azimuth_phi_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_phi[_j][_qp](2))*std::sin(2.*_polar_theta[_qp])))/((1. + Utility::pow<2>(_alpha))*_M));
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
      return _test[_i][_qp] * ((-0.5*_g0*_phi[_j][_qp]*Utility::pow<3>((1.0/std::sin(_polar_theta[_qp])))*(-4.*_alpha*_Ku*Utility::pow<2>(_M)*_nx*_ny*std::sin(2.*_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]) + 2.*_Ku*Utility::pow<2>(_M)*_ny*_nz*std::sin(2.*_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]) + 4.*_Ae*Utility::pow<2>(_azimuth_phi_grad[_qp](0))*Utility::pow<4>(std::sin(_polar_theta[_qp])) + 
       4.*_Ae*Utility::pow<2>(_azimuth_phi_grad[_qp](1))*Utility::pow<4>(std::sin(_polar_theta[_qp])) + 4.*_Ae*Utility::pow<2>(_azimuth_phi_grad[_qp](2))*Utility::pow<4>(std::sin(_polar_theta[_qp])) - 4.*_Ku*Utility::pow<2>(_M)*Utility::pow<2>(_nx)*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*Utility::pow<4>(std::sin(_polar_theta[_qp])) + 
       4.*_Ku*Utility::pow<2>(_M)*Utility::pow<2>(_ny)*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*Utility::pow<4>(std::sin(_polar_theta[_qp])) + 4.*_Ku*Utility::pow<2>(_M)*std::cos(_polar_theta[_qp])*(_alpha*(-1.*Utility::pow<2>(_nx) + Utility::pow<2>(_nz))*std::sin(2.*_azimuth_phi[_qp]) + 4.*_nx*_ny*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp]))) - 
       2.*_alpha*_Ku*Utility::pow<2>(_M)*_nz*std::cos(2.*_azimuth_phi[_qp])*(_nx*(3. + std::cos(2.*_polar_theta[_qp])) + _ny*std::sin(2.*_polar_theta[_qp])) + _Ku*Utility::pow<2>(_M)*_nx*_ny*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*std::sin(4.*_polar_theta[_qp])))/((1. + Utility::pow<2>(_alpha))*_M));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
