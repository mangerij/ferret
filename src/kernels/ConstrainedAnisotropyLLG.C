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

#include "ConstrainedAnisotropyLLG.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<ConstrainedAnisotropyLLG>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for theta, 1 for phi)");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  params.addRequiredParam<Real>("alpha", "the damping coefficient in the LLG equation");
  params.addRequiredParam<Real>("Ku", "Ku");
  params.addRequiredParam<Real>("nx", "nx");
  params.addRequiredParam<Real>("ny", "ny");
  params.addRequiredParam<Real>("nz", "nz");
  params.addRequiredParam<Real>("g0", "g0");
  params.addRequiredParam<Real>("M", "M");
  return params;
}

ConstrainedAnisotropyLLG::ConstrainedAnisotropyLLG(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _azimuth_phi_var(coupled("azimuth_phi")),
  _polar_theta_var(coupled("polar_theta")),
  _azimuth_phi(coupledValue("azimuth_phi")),
  _polar_theta(coupledValue("polar_theta")),
  _azimuth_phi_grad(coupledGradient("azimuth_phi")),
  _polar_theta_grad(coupledGradient("polar_theta")),
  _alpha(getParam<Real>("alpha")),
  _Ku(getParam<Real>("Ku")),
  _nx(getParam<Real>("nx")),
  _ny(getParam<Real>("ny")),
  _nz(getParam<Real>("nz")),
  _g0(getParam<Real>("g0")),
  _M(getParam<Real>("M"))
{
}

Real
ConstrainedAnisotropyLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * ((-2*_g0*_Ku*_M*(std::cos(_azimuth_phi[_qp])*(_ny + _alpha*_nx*std::cos(_polar_theta[_qp])) + (-_nx + _alpha*_ny*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) - _alpha*_nz*std::sin(_polar_theta[_qp]))*(_nz*std::cos(_polar_theta[_qp]) + (_nx*std::cos(_azimuth_phi[_qp]) + _ny*std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp])))/(1.0 + Utility::pow<2>(_alpha)));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * ((2*_g0*_Ku*_M*(_nx*std::cos(_azimuth_phi[_qp]) + _nz*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])) + _ny*std::sin(_azimuth_phi[_qp]))*(-_nz + std::cos(_azimuth_phi[_qp])*(-(_alpha*_ny) + _nx*std::cos(_polar_theta[_qp]))*(1.0/std::sin(_polar_theta[_qp])) + _ny*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + _alpha*_nx*(1.0/std::sin(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp]))/(1.0 + Utility::pow<2>(_alpha)));
  }
  else
    return 0.0;
}

Real
ConstrainedAnisotropyLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp]*((_g0*_Ku*_M*_phi[_j][_qp]*(-(_alpha*std::cos(2*_polar_theta[_qp])*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny) - 2*Utility::pow<2>(_nz) + (_nx - _ny)*(_nx + _ny)*std::cos(2*_azimuth_phi[_qp]) + 2*_nx*_ny*std::sin(2*_azimuth_phi[_qp]))) + std::cos(_polar_theta[_qp])*(-2*_nx*_ny*std::cos(2*_azimuth_phi[_qp]) + (_nx - _ny)*(_nx + _ny)*std::sin(2*_azimuth_phi[_qp])) + 2*_nz*(_ny*std::cos(_azimuth_phi[_qp]) - _nx*std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp]) + 4*_alpha*_nz*(_nx*std::cos(_azimuth_phi[_qp]) + _ny*std::sin(_azimuth_phi[_qp]))*std::sin(2*_polar_theta[_qp])))/(1 + Utility::pow<2>(_alpha)));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp]*((2*_g0*_Ku*_M*_phi[_j][_qp]*(Utility::pow<2>(std::cos(_azimuth_phi[_qp]))*(_alpha*(_nx - _ny)*(_nx + _ny) + 2*_nx*_ny*std::cos(_polar_theta[_qp])) + 
       std::sin(_azimuth_phi[_qp])*(_nz*(_alpha*_ny - _nx*std::cos(_polar_theta[_qp]))*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])) + (-(_alpha*Utility::pow<2>(_nx)) + _alpha*Utility::pow<2>(_ny) - 2*_nx*_ny*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + _nx*_nz*std::sin(_polar_theta[_qp])) + 
       std::cos(_azimuth_phi[_qp])*(_nz*(_alpha*_nx + _ny*std::cos(_polar_theta[_qp]))*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])) + 2*(2*_alpha*_nx*_ny + (-Utility::pow<2>(_nx) + Utility::pow<2>(_ny))*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) - _ny*_nz*std::sin(_polar_theta[_qp]))))/(1.0 + Utility::pow<2>(_alpha)));
  }
  else
    return 0.0;
}

Real
ConstrainedAnisotropyLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuth_phi_var)
    {
      return _test[_i][_qp]*((2*_g0*_Ku*_M*_phi[_j][_qp]*(_nz*std::cos(_azimuth_phi[_qp])*(_nx*std::cos(_polar_theta[_qp]) - _alpha*_ny*std::cos(2*_polar_theta[_qp])) + _nz*(_ny*std::cos(_polar_theta[_qp]) + _alpha*_nx*std::cos(2*_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + std::cos(2*_azimuth_phi[_qp])*(Utility::pow<2>(_nx) - Utility::pow<2>(_ny) - 2*_alpha*_nx*_ny*std::cos(_polar_theta[_qp]))*std::sin(_polar_theta[_qp]) + 
       (2*_nx*_ny + _alpha*(_nx - _ny)*(_nx + _ny)*std::cos(_polar_theta[_qp]))*std::sin(2*_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])))/(1.0 + Utility::pow<2>(_alpha)));
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
      return _test[_i][_qp] * (-((_g0*_Ku*_M*_phi[_j][_qp]*(-2*Utility::pow<2>(_nz) + 2*Utility::pow<2>(_nx)*Utility::pow<2>(std::cos(_azimuth_phi[_qp])) - _nz*std::cos(_azimuth_phi[_qp])*(2*_alpha*_ny - 3*_nx*std::cos(_polar_theta[_qp]) + _nx*std::cos(3*_polar_theta[_qp]))*Utility::pow<3>((1.0/std::sin(_polar_theta[_qp]))) + 
         2*std::sin(_azimuth_phi[_qp])*(_alpha*_nx*_nz*Utility::pow<3>((1.0/std::sin(_polar_theta[_qp]))) + _ny*_nz*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*(2 + Utility::pow<2>((1.0/std::sin(_polar_theta[_qp])))) + Utility::pow<2>(_ny)*std::sin(_azimuth_phi[_qp])) + 2*_nx*_ny*std::sin(2*_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp]))/(1 + Utility::pow<2>(_alpha))));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
