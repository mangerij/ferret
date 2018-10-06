/**
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT any WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "AnisotropyConstrainedLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", AnisotropyConstrainedLLG);

template<>
InputParameters validParams<AnisotropyConstrainedLLG>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for theta, 1 for phi)");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  params.addRequiredParam<Real>("alpha", "the damping coefficient in the LLG equation");
  params.addRequiredParam<Real>("K1", "K1");
  params.addRequiredParam<Real>("K2", "K2");
  params.addRequiredParam<Real>("nx", "nx");
  params.addRequiredParam<Real>("ny", "ny");
  params.addRequiredParam<Real>("nz", "nz");
  params.addRequiredParam<Real>("g0", "g0");
  params.addRequiredParam<Real>("M", "M");
  return params;
}

AnisotropyConstrainedLLG::AnisotropyConstrainedLLG(const InputParameters & parameters)
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
  _nx(getParam<Real>("nx")),
  _ny(getParam<Real>("ny")),
  _nz(getParam<Real>("nz")),
  _g0(getParam<Real>("g0"))
{
}

Real
AnisotropyConstrainedLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return (-(_g0*_test[_i][_qp]*((1 + Utility::pow<2>(_alpha))*(_ny*std::cos(_azimuth_phi[_qp]) - _nx*std::sin(_azimuth_phi[_qp])) + _alpha*std::cos(_polar_theta[_qp])*(_nx*std::cos(_azimuth_phi[_qp]) + _ny*std::sin(_azimuth_phi[_qp])) - _alpha*_nz*std::sin(_polar_theta[_qp]))*
      (((8*_K1*(_nx + _ny) + 3*_K2*_nx*(3*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) + 4*Utility::pow<2>(_nz)))*std::cos(_azimuth_phi[_qp]) + 
           3*_K2*(_nx*(Utility::pow<2>(_nx) - 3*Utility::pow<2>(_ny))*std::cos(3*_azimuth_phi[_qp]) + _ny*(3*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) + 4*Utility::pow<2>(_nz))*std::sin(_azimuth_phi[_qp]) - _ny*(-3*Utility::pow<2>(_nx) + Utility::pow<2>(_ny))*std::sin(3*_azimuth_phi[_qp])))*std::sin(_polar_theta[_qp]) + 
        4*_nz*std::cos(_polar_theta[_qp])*(2*_K1 + 3*_K2*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) + 2*_K2*Utility::pow<2>(_nz) + _K2*(-3*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) + 2*Utility::pow<2>(_nz))*std::cos(2*_polar_theta[_qp]) + 
           6*_K2*((_nx - _ny)*(_nx + _ny)*std::cos(2*_azimuth_phi[_qp]) + 2*_nx*_ny*std::sin(2*_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp]))) - 
        2*_K2*(_nx*std::cos(_azimuth_phi[_qp]) + _ny*std::sin(_azimuth_phi[_qp]))*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny) - 6*Utility::pow<2>(_nz) + (_nx - _ny)*(_nx + _ny)*std::cos(2*_azimuth_phi[_qp]) + 2*_nx*_ny*std::sin(2*_azimuth_phi[_qp]))*std::sin(3*_polar_theta[_qp])))/(4.*(1 + Utility::pow<2>(_alpha))));
  }
  else if (_component == 1)
  {
    return ((2*_g0*_test[_i][_qp]*(-((1 + Utility::pow<2>(_alpha))*_nz) + std::cos(_azimuth_phi[_qp])*((1 + Utility::pow<2>(_alpha))*_nx*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])) - _alpha*_ny*(1.0/std::sin(_polar_theta[_qp]))) + (1 + Utility::pow<2>(_alpha))*_ny*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + _alpha*_nx*(1.0/std::sin(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]))*
     (_nz*(_K1 + 2*_K2*Utility::pow<2>(_nz))*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*Utility::pow<2>((1.0/std::sin(_polar_theta[_qp]))) + Utility::pow<2>((1.0/std::sin(_polar_theta[_qp])))*((_K1*(_nx + _ny) + 6*_K2*_nx*Utility::pow<2>(_nz))*std::cos(_azimuth_phi[_qp]) + 6*_K2*_ny*Utility::pow<2>(_nz)*std::sin(_azimuth_phi[_qp])) + 
       _K2*(_nx*std::cos(_azimuth_phi[_qp]) + _ny*std::sin(_azimuth_phi[_qp]))*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny) - 6*Utility::pow<2>(_nz) + (_nx - _ny)*(_nx + _ny)*std::cos(2*_azimuth_phi[_qp]) + 2*_nx*_ny*std::sin(2*_azimuth_phi[_qp])) + 
       _K2*_nz*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*(3*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) - 2*Utility::pow<2>(_nz) + 3*(_nx - _ny)*(_nx + _ny)*std::cos(2*_azimuth_phi[_qp]) + 6*_nx*_ny*std::sin(2*_azimuth_phi[_qp])))*Utility::pow<3>(std::sin(_polar_theta[_qp])))/(1 + Utility::pow<2>(_alpha)));
  }
  else
    return 0.0;
}

Real
AnisotropyConstrainedLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return ((_g0*_phi[_j][_qp]*_test[_i][_qp]*(3*(1 + Utility::pow<2>(_alpha))*_K2*std::cos(3*_polar_theta[_qp])*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny) - 6*Utility::pow<2>(_nz) + (_nx - _ny)*(_nx + _ny)*std::cos(2*_azimuth_phi[_qp]) + 2*_nx*_ny*std::sin(2*_azimuth_phi[_qp]))*(2*_nx*_ny*std::cos(2*_azimuth_phi[_qp]) + (-Utility::pow<2>(_nx) + Utility::pow<2>(_ny))*std::sin(2*_azimuth_phi[_qp])) + 
       _alpha*_K2*std::cos(4*_polar_theta[_qp])*(3*Utility::pow<2>(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) - 24*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny))*Utility::pow<2>(_nz) + 8*Utility::pow<4>(_nz) + (Utility::pow<4>(_nx) - 6*Utility::pow<2>(_nx)*Utility::pow<2>(_ny) + Utility::pow<4>(_ny))*std::cos(4*_azimuth_phi[_qp]) + 
          8*_nx*_ny*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny) - 6*Utility::pow<2>(_nz))*std::sin(2*_azimuth_phi[_qp]) + 4*(_nx - _ny)*(_nx + _ny)*std::cos(2*_azimuth_phi[_qp])*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny) - 6*Utility::pow<2>(_nz) + 2*_nx*_ny*std::sin(2*_azimuth_phi[_qp]))) - 
       (1 + Utility::pow<2>(_alpha))*std::cos(_polar_theta[_qp])*(_ny*std::cos(_azimuth_phi[_qp]) - _nx*std::sin(_azimuth_phi[_qp]))*((8*_K1*(_nx + _ny) + 3*_K2*_nx*(3*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) + 4*Utility::pow<2>(_nz)))*std::cos(_azimuth_phi[_qp]) + 
          3*_K2*(_nx*(Utility::pow<2>(_nx) - 3*Utility::pow<2>(_ny))*std::cos(3*_azimuth_phi[_qp]) + _ny*(3*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) + 4*Utility::pow<2>(_nz))*std::sin(_azimuth_phi[_qp]) - _ny*(-3*Utility::pow<2>(_nx) + Utility::pow<2>(_ny))*std::sin(3*_azimuth_phi[_qp]))) - 
       _alpha*std::cos(2*_polar_theta[_qp])*(4*_K1*_nx*(_nx + _ny) + 3*_K2*Utility::pow<2>(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) - 8*_K1*Utility::pow<2>(_nz) - 8*_K2*Utility::pow<4>(_nz) + 4*(_K1*_nx*(_nx + _ny) + _K2*(Utility::pow<4>(_nx) - Utility::pow<4>(_ny)))*std::cos(2*_azimuth_phi[_qp]) + 
          _K2*(Utility::pow<4>(_nx) - 6*Utility::pow<2>(_nx)*Utility::pow<2>(_ny) + Utility::pow<4>(_ny))*std::cos(4*_azimuth_phi[_qp]) + 4*_ny*(_K1*(_nx + _ny) + 2*_K2*_nx*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)))*std::sin(2*_azimuth_phi[_qp]) + 4*_K2*_nx*(_nx - _ny)*_ny*(_nx + _ny)*std::sin(4*_azimuth_phi[_qp])) + 
       2*(1 + Utility::pow<2>(_alpha))*_nz*(_ny*std::cos(_azimuth_phi[_qp]) - _nx*std::sin(_azimuth_phi[_qp]))*(4*_K1 + 3*_K2*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny) + 2*Utility::pow<2>(_nz)) + 3*_K2*(_nx - _ny)*(_nx + _ny)*std::cos(2*_azimuth_phi[_qp]) + 6*_K2*_nx*_ny*std::sin(2*_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp]) + 
       4*_alpha*_nz*((2*_K1*(2*_nx + _ny) + _K2*_nx*(3*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) + 4*Utility::pow<2>(_nz)))*std::cos(_azimuth_phi[_qp]) + _K2*_nx*(Utility::pow<2>(_nx) - 3*Utility::pow<2>(_ny))*std::cos(3*_azimuth_phi[_qp]) + 
          _ny*(2*_K1 + 3*_K2*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) + 4*_K2*Utility::pow<2>(_nz))*std::sin(_azimuth_phi[_qp]) - _K2*_ny*(-3*Utility::pow<2>(_nx) + Utility::pow<2>(_ny))*std::sin(3*_azimuth_phi[_qp]))*std::sin(2*_polar_theta[_qp]) - 
       6*(1 + Utility::pow<2>(_alpha))*_K2*_nz*(_ny*std::cos(_azimuth_phi[_qp]) - _nx*std::sin(_azimuth_phi[_qp]))*(3*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) - 2*Utility::pow<2>(_nz) + 3*(_nx - _ny)*(_nx + _ny)*std::cos(2*_azimuth_phi[_qp]) + 6*_nx*_ny*std::sin(2*_azimuth_phi[_qp]))*std::sin(3*_polar_theta[_qp]) - 
       16*_alpha*_K2*_nz*(_nx*std::cos(_azimuth_phi[_qp]) + _ny*std::sin(_azimuth_phi[_qp]))*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny) - 2*Utility::pow<2>(_nz) + (_nx - _ny)*(_nx + _ny)*std::cos(2*_azimuth_phi[_qp]) + 2*_nx*_ny*std::sin(2*_azimuth_phi[_qp]))*std::sin(4*_polar_theta[_qp])))/(4.*(1 + Utility::pow<2>(_alpha))));
  }
  else if (_component == 1)
  {
    return ((2*_g0*_phi[_j][_qp]*_test[_i][_qp]*(1.0/std::sin(_polar_theta[_qp]))*(-2*(1 + Utility::pow<2>(_alpha))*_K2*_nx*Utility::pow<3>(_nz)*Utility::pow<4>(std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) - _alpha*_K1*Utility::pow<2>(_nx)*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp]) - _alpha*_K1*_nx*_ny*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp]) + 
       _alpha*_K1*_nx*_ny*std::sin(2*_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]) + _alpha*_K1*Utility::pow<2>(_ny)*std::sin(2*_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]) + _K1*_nx*_nz*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) + Utility::pow<2>(_alpha)*_K1*_nx*_nz*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 
       _K1*_ny*_nz*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) + Utility::pow<2>(_alpha)*_K1*_ny*_nz*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) - 6*_alpha*_K2*Utility::pow<2>(_nx)*Utility::pow<2>(_ny)*Utility::pow<4>(std::sin(_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp])) + 
       2*_alpha*_K2*Utility::pow<4>(_ny)*Utility::pow<4>(std::sin(_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp])) + 9*_alpha*_K2*Utility::pow<2>(_nx)*Utility::pow<2>(_ny)*Utility::pow<2>(std::sin(2*_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp])) - 
       6*(1 + Utility::pow<2>(_alpha))*_K2*Utility::pow<2>(_nx)*_ny*_nz*Utility::pow<5>(std::cos(_azimuth_phi[_qp]))*Utility::pow<4>(std::sin(_polar_theta[_qp])) + 6*_K2*_nx*Utility::pow<2>(_ny)*_nz*Utility::pow<5>(std::sin(_azimuth_phi[_qp]))*Utility::pow<4>(std::sin(_polar_theta[_qp])) + 
       6*Utility::pow<2>(_alpha)*_K2*_nx*Utility::pow<2>(_ny)*_nz*Utility::pow<5>(std::sin(_azimuth_phi[_qp]))*Utility::pow<4>(std::sin(_polar_theta[_qp])) + 2*_K2*Utility::pow<2>(_nz)*Utility::pow<3>(std::cos(_polar_theta[_qp]))*
        (_alpha*_ny*_nz*std::sin(_azimuth_phi[_qp]) - 3*(1 + Utility::pow<2>(_alpha))*(2*_nx*_ny*Utility::pow<2>(std::sin(_azimuth_phi[_qp])) + (_nx - _ny)*(_nx + _ny)*std::sin(2*_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp])) + 
       2*_K2*Utility::pow<3>(std::cos(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp]))*(3*_nz*std::cos(_polar_theta[_qp])*(_alpha*_nx*(Utility::pow<2>(_nx) - 2*Utility::pow<2>(_ny)) + (1 + Utility::pow<2>(_alpha))*_ny*(3*Utility::pow<2>(_nx) - Utility::pow<2>(_nz))*std::cos(_polar_theta[_qp])) + 
          2*(_alpha*_nx*_ny*(5*Utility::pow<2>(_nx) - 3*Utility::pow<2>(_ny)) - (1 + Utility::pow<2>(_alpha))*(2*Utility::pow<4>(_nx) + 3*Utility::pow<2>(_ny)*Utility::pow<2>(_nz) - 3*Utility::pow<2>(_nx)*(2*Utility::pow<2>(_ny) + Utility::pow<2>(_nz)))*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]) + 
          3*(1 + Utility::pow<2>(_alpha))*(_nx - _ny)*_ny*(_nx + _ny)*_nz*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp]))) + 
       std::cos(_polar_theta[_qp])*(_alpha*_K1*_ny*_nz*std::sin(_azimuth_phi[_qp]) - (1 + Utility::pow<2>(_alpha))*_K1*_ny*(_nx + _ny)*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp]) + 
          6*_alpha*_K2*_nz*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*(-2*Utility::pow<3>(_nx)*std::cos(_azimuth_phi[_qp]) + _ny*(-2*Utility::pow<2>(_nx) + Utility::pow<2>(_ny))*std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) - 
          2*(1 + Utility::pow<2>(_alpha))*_K2*_nx*_ny*((4*Utility::pow<2>(_ny) - 6*Utility::pow<2>(_nz))*Utility::pow<4>(std::sin(_azimuth_phi[_qp])) + 3*(_nx - _ny)*(_nx + _ny)*Utility::pow<2>(std::sin(2*_azimuth_phi[_qp])))*Utility::pow<3>(std::sin(_polar_theta[_qp]))) + 
       _nz*Utility::pow<2>(std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])*(-((1 + Utility::pow<2>(_alpha))*_K1*_nx) + 6*_K2*_nz*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])*(_alpha*(-Utility::pow<2>(_nx) + Utility::pow<2>(_ny)) + (1 + Utility::pow<2>(_alpha))*_nx*_nz*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]))) + 
       2*_K2*_nx*Utility::pow<4>(std::cos(_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp]))*(_alpha*_nx*(Utility::pow<2>(_nx) - 3*Utility::pow<2>(_ny)) + 
          (1 + Utility::pow<2>(_alpha))*(2*_ny*(2*Utility::pow<2>(_nx) - 3*Utility::pow<2>(_nz))*std::cos(_polar_theta[_qp]) + 3*(Utility::pow<2>(_nx) - 2*Utility::pow<2>(_ny))*_nz*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]))) + 
       Utility::pow<2>(std::cos(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp])*(_alpha*_K1*_nx*(_nx + _ny) + 12*(1 + Utility::pow<2>(_alpha))*_K2*_nx*_ny*Utility::pow<2>(_nz)*Utility::pow<3>(std::cos(_polar_theta[_qp])) + 
          _ny*std::cos(_polar_theta[_qp])*((1 + Utility::pow<2>(_alpha))*_K1*(_nx + _ny) + 6*_alpha*_K2*(7*Utility::pow<2>(_nx) - 2*Utility::pow<2>(_ny))*_nz*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])) + 
          6*_K2*_nz*Utility::pow<2>(std::cos(_polar_theta[_qp]))*(_alpha*(_nx - _ny)*(_nx + _ny)*_nz + (1 + Utility::pow<2>(_alpha))*_nx*(-3*Utility::pow<2>(_nx) + 6*Utility::pow<2>(_ny) + Utility::pow<2>(_nz))*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])) + 
          std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])*((1 + Utility::pow<2>(_alpha))*_K1*(_nx + _ny)*_nz + 6*_K2*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])*(-(_alpha*(Utility::pow<4>(_nx) + Utility::pow<4>(_ny))) + (1 + Utility::pow<2>(_alpha))*_nx*(_nx - _ny)*(_nx + _ny)*_nz*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])))) - 
       (9*(1 + Utility::pow<2>(_alpha))*_K2*_nx*Utility::pow<2>(_ny)*_nz*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(2*_polar_theta[_qp])))/2. + 
       (std::cos(_azimuth_phi[_qp])*(4*_alpha*_K2*_nx*Utility::pow<3>(_nz)*Utility::pow<3>(std::cos(_polar_theta[_qp])) + 4*(1 + Utility::pow<2>(_alpha))*_K2*_ny*Utility::pow<3>(_nz)*Utility::pow<4>(std::cos(_polar_theta[_qp])) + 2*_ny*_nz*Utility::pow<2>(std::cos(_polar_theta[_qp]))*(_K1 + Utility::pow<2>(_alpha)*_K1 + 24*_alpha*_K2*_nx*_nz*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])) + 
            2*std::cos(_polar_theta[_qp])*(_alpha*_K1*_nx*_nz + 2*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])*(-((1 + Utility::pow<2>(_alpha))*_K1*_nx*(_nx + _ny)) + 
                  _K2*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])*(21*_alpha*_nx*Utility::pow<2>(_ny)*_nz + 2*(1 + Utility::pow<2>(_alpha))*(2*Utility::pow<4>(_ny) - 3*Utility::pow<2>(_ny)*Utility::pow<2>(_nz) + 3*Utility::pow<2>(_nx)*(-2*Utility::pow<2>(_ny) + Utility::pow<2>(_nz)))*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])))) + 
            _K2*_ny*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*(4*std::sin(_azimuth_phi[_qp])*Utility::pow<3>(std::sin(_polar_theta[_qp]))*(-6*_alpha*Utility::pow<3>(_nx) + 10*_alpha*_nx*Utility::pow<2>(_ny) + 3*(1 + Utility::pow<2>(_alpha))*(2*Utility::pow<2>(_nx) - Utility::pow<2>(_ny))*_nz*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])) - 
               3*(1 + Utility::pow<2>(_alpha))*_nz*(6*Utility::pow<2>(_nx) - 3*Utility::pow<2>(_ny) + Utility::pow<2>(_nz))*Utility::pow<2>(std::sin(2*_polar_theta[_qp])))))/2.))/(1 + Utility::pow<2>(_alpha)));
  }
  else
    return 0.0;
}

Real
AnisotropyConstrainedLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuth_phi_var)
    {
      return ((2*_g0*_phi[_j][_qp]*_test[_i][_qp]*(6*_K2*_ny*_nz*(-((1 + Utility::pow<2>(_alpha))*(2*Utility::pow<2>(_nx) - Utility::pow<2>(_ny))*std::cos(_polar_theta[_qp])) + _alpha*_nx*_ny*(1 + 2*std::cos(2*_polar_theta[_qp])))*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 
       2*_K2*Utility::pow<2>(_nx)*Utility::pow<4>(std::cos(_azimuth_phi[_qp]))*((1 + Utility::pow<2>(_alpha))*(Utility::pow<2>(_nx) - 3*Utility::pow<2>(_ny)) - 4*_alpha*_nx*_ny*std::cos(_polar_theta[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp])) + 
       2*_K2*Utility::pow<2>(_ny)*(-((1 + Utility::pow<2>(_alpha))*(3*Utility::pow<2>(_nx) - Utility::pow<2>(_ny))) + 4*_alpha*_nx*_ny*std::cos(_polar_theta[_qp]))*Utility::pow<4>(std::sin(_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp])) + 
       2*_K2*_nx*Utility::pow<3>(std::cos(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp]))*(3*_alpha*_nx*_ny*_nz + 3*_nz*std::cos(_polar_theta[_qp])*((1 + Utility::pow<2>(_alpha))*(Utility::pow<2>(_nx) - 2*Utility::pow<2>(_ny)) - 4*_alpha*_nx*_ny*std::cos(_polar_theta[_qp])) + 
          2*((1 + Utility::pow<2>(_alpha))*_ny*(5*Utility::pow<2>(_nx) - 3*Utility::pow<2>(_ny)) + 2*_alpha*_nx*(Utility::pow<2>(_nx) - 3*Utility::pow<2>(_ny))*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])) + 
       _nz*std::sin(_azimuth_phi[_qp])*(-(_alpha*_K1*(_nx + _ny)) + std::cos(_polar_theta[_qp])*(_alpha*_K1*(2*_nx + _ny)*std::cos(_polar_theta[_qp]) + (1 + Utility::pow<2>(_alpha))*_ny*(_K1 + _K2*Utility::pow<2>(_nz) + _K2*Utility::pow<2>(_nz)*std::cos(2*_polar_theta[_qp])) + 2*_alpha*_K2*_nx*Utility::pow<2>(_nz)*std::cos(3*_polar_theta[_qp])) + 
          6*_K2*Utility::pow<2>(_nx)*(-((1 + Utility::pow<2>(_alpha))*_nx*std::cos(_polar_theta[_qp])) + _alpha*_ny*(1 + 2*std::cos(2*_polar_theta[_qp])))*std::sin(2*_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp]))) - 
       Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*(((1 + Utility::pow<2>(_alpha))*_K1*_nx*(_nx + _ny) + std::cos(_polar_theta[_qp])*(-(_alpha*_K1*_ny*(_nx + _ny)) + 6*_K2*Utility::pow<2>(_nz)*std::cos(_polar_theta[_qp])*((1 + Utility::pow<2>(_alpha))*(_nx - _ny)*(_nx + _ny) - 4*_alpha*_nx*_ny*std::cos(_polar_theta[_qp]))))*std::sin(_polar_theta[_qp]) + 
          6*_alpha*_K2*_nx*_ny*Utility::pow<2>(_nz)*std::sin(2*_polar_theta[_qp])) + std::sin(2*_azimuth_phi[_qp])*((1 + Utility::pow<2>(_alpha))*_K1*_ny*(_nx + _ny)*std::sin(_polar_theta[_qp]) + 
          3*_K2*_nx*_ny*(3*(1 + Utility::pow<2>(_alpha))*_nx*_ny + 2*_alpha*(_nx - _ny)*(_nx + _ny)*std::cos(_polar_theta[_qp]))*std::sin(2*_azimuth_phi[_qp])*Utility::pow<3>(std::sin(_polar_theta[_qp])) + 3*_alpha*_K2*(-Utility::pow<2>(_nx) + Utility::pow<2>(_ny))*Utility::pow<2>(_nz)*std::sin(2*_polar_theta[_qp])) + 
       std::cos(_azimuth_phi[_qp])*(-8*_alpha*_K2*_ny*Utility::pow<3>(_nz)*Utility::pow<4>(std::cos(_polar_theta[_qp])) + 2*_K2*Utility::pow<2>(_nz)*Utility::pow<3>(std::cos(_polar_theta[_qp]))*((1 + Utility::pow<2>(_alpha))*_nx*_nz + 12*_alpha*(_nx - _ny)*(_nx + _ny)*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])) + 
          _ny*_nz*Utility::pow<2>(std::cos(_polar_theta[_qp]))*(-(_alpha*(_K1 - 6*_K2*Utility::pow<2>(_nz))) + 24*(1 + Utility::pow<2>(_alpha))*_K2*_nx*_nz*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])) + 
          std::cos(_polar_theta[_qp])*((1 + Utility::pow<2>(_alpha))*_K1*_nx*_nz + 2*_K2*Utility::pow<2>(_ny)*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp]))*(21*(1 + Utility::pow<2>(_alpha))*_nx*_nz - 4*_alpha*(-3*Utility::pow<2>(_nx) + Utility::pow<2>(_ny))*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]))) + 
          std::sin(_azimuth_phi[_qp])*(2*_K2*std::sin(_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp]))*(-3*_alpha*Utility::pow<3>(_ny)*_nz*(1 + 2*std::cos(2*_polar_theta[_qp])) + 2*(1 + Utility::pow<2>(_alpha))*_nx*(-3*Utility::pow<2>(_nx)*_ny + 5*Utility::pow<3>(_ny))*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])) + _alpha*_K1*_nx*(_nx + _ny)*std::sin(2*_polar_theta[_qp]))) + 
       Utility::pow<2>(std::cos(_azimuth_phi[_qp]))*(6*_K2*_nx*_nz*(-(_alpha*(Utility::pow<2>(_nx) - 2*Utility::pow<2>(_ny))) + 7*(1 + Utility::pow<2>(_alpha))*_nx*_ny*std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp])) - 
          6*(1 + Utility::pow<2>(_alpha))*_K2*(Utility::pow<4>(_nx) + Utility::pow<4>(_ny))*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp])) + 6*_alpha*_K2*_nx*_nz*std::sin(2*_polar_theta[_qp])*(_ny*_nz + (Utility::pow<2>(_nx) - 2*Utility::pow<2>(_ny))*std::sin(_azimuth_phi[_qp])*std::sin(2*_polar_theta[_qp])) + 
          std::sin(_polar_theta[_qp])*((1 + Utility::pow<2>(_alpha))*_K1*_nx*(_nx + _ny) + std::cos(_polar_theta[_qp])*(-(_alpha*_K1*_ny*(_nx + _ny)) + 6*_K2*Utility::pow<2>(_nz)*std::cos(_polar_theta[_qp])*((1 + Utility::pow<2>(_alpha))*(_nx - _ny)*(_nx + _ny) - 4*_alpha*_nx*_ny*std::cos(_polar_theta[_qp]))) - 
             6*(1 + Utility::pow<2>(_alpha))*_K2*Utility::pow<3>(_ny)*_nz*std::sin(_azimuth_phi[_qp])*std::sin(2*_polar_theta[_qp])))))/(1 + Utility::pow<2>(_alpha)));
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
      return ((-2*_g0*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>((1.0/std::sin(_polar_theta[_qp])))*(2*_alpha*_K2*_nx*Utility::pow<3>(_nz)*Utility::pow<4>(std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + 2*(1 + Utility::pow<2>(_alpha))*_K2*_ny*Utility::pow<3>(_nz)*Utility::pow<5>(std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + 
       6*(1 + Utility::pow<2>(_alpha))*_K2*Utility::pow<3>(_nx)*_nz*Utility::pow<5>(std::cos(_azimuth_phi[_qp]))*std::cos(_polar_theta[_qp])*Utility::pow<4>(std::sin(_polar_theta[_qp])) + 
       2*_K2*std::cos(_polar_theta[_qp])*Utility::pow<3>(std::sin(_polar_theta[_qp]))*(-2*_alpha*_nx*_ny*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*(-3*Utility::pow<2>(_nz) + Utility::pow<2>(_ny)*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))) + 3*_alpha*(_nx - _ny)*(_nx + _ny)*Utility::pow<2>(_nz)*std::sin(2*_azimuth_phi[_qp]) + 
          3*(1 + Utility::pow<2>(_alpha))*_ny*_nz*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*(2*(_ny - _nz)*(_ny + _nz) + Utility::pow<2>(_ny)*Utility::pow<2>(std::sin(_azimuth_phi[_qp])))*std::sin(_polar_theta[_qp])) + 
       Utility::pow<3>(std::cos(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp]))*((_nz*((1 + Utility::pow<2>(_alpha))*(4*_K1*(_nx + _ny) - 3*_K2*_nx*(Utility::pow<2>(_nx) - 3*Utility::pow<2>(_ny) - 2*Utility::pow<2>(_nz)))*std::cos(_polar_theta[_qp]) + 
               3*_K2*_nx*(8*_alpha*_nx*_ny*std::cos(2*_polar_theta[_qp]) - (1 + Utility::pow<2>(_alpha))*(7*Utility::pow<2>(_nx) + 3*Utility::pow<2>(_ny) - 6*Utility::pow<2>(_nz))*std::cos(3*_polar_theta[_qp]))))/4. - 
          2*_K2*_nx*(2*_alpha*_nx*(Utility::pow<2>(_nx) - 3*Utility::pow<2>(_ny))*std::cos(_polar_theta[_qp]) + (1 + Utility::pow<2>(_alpha))*_ny*(2*Utility::pow<2>(_nx) - 3*Utility::pow<2>(_nz))*(1 + 3*std::cos(2*_polar_theta[_qp])))*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]) - 
          3*(1 + Utility::pow<2>(_alpha))*_K2*_nx*(Utility::pow<2>(_nx) + 3*Utility::pow<2>(_ny))*_nz*std::cos(2*_azimuth_phi[_qp])*std::cos(_polar_theta[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp]))) + 
       (1 + Utility::pow<2>(_alpha))*_ny*_nz*Utility::pow<3>(std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])*(_K1 + 6*_K2*(-Utility::pow<2>(_ny) + Utility::pow<2>(_nz))*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp]))) + 
       (Utility::pow<2>(std::cos(_polar_theta[_qp]))*(_alpha*_nx*_nz*(4*_K1 - 9*_K2*Utility::pow<2>(_ny) + 12*_K2*Utility::pow<2>(_nz) + 3*_K2*(3*Utility::pow<2>(_ny) - 4*Utility::pow<2>(_nz))*std::cos(2*_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) - 
            2*(1 + Utility::pow<2>(_alpha))*_K2*(-2*(2*Utility::pow<4>(_ny) - 9*Utility::pow<2>(_ny)*Utility::pow<2>(_nz) + 3*Utility::pow<4>(_nz))*std::cos(2*_azimuth_phi[_qp]) + 
               (Utility::pow<4>(_ny) + 3*Utility::pow<2>(_ny)*Utility::pow<2>(_nz) + 3*Utility::pow<2>(_nx)*(-2*Utility::pow<2>(_ny) + Utility::pow<2>(_nz)))*std::cos(4*_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp]))))/4. + 
       (32*(1 + Utility::pow<2>(_alpha))*_ny*(_K1*(_nx + _ny) + _K2*_nx*(-2*Utility::pow<2>(_ny) + 21*Utility::pow<2>(_nz)) + 3*_K2*_nx*(-2*Utility::pow<2>(_ny) + 9*Utility::pow<2>(_nz))*std::cos(2*_polar_theta[_qp]) + _K2*_nx*(2*Utility::pow<2>(_ny) - 3*Utility::pow<2>(_nz))*std::cos(2*_azimuth_phi[_qp])*(1 + 3*std::cos(2*_polar_theta[_qp])))*
           std::sin(2*_azimuth_phi[_qp])*Utility::pow<3>(std::sin(_polar_theta[_qp])) + std::cos(_azimuth_phi[_qp])*(-16*_alpha*_ny*_nz*(4*_K1 + 3*_K2*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny) + 2*Utility::pow<2>(_nz))) + 
             2*std::cos(_polar_theta[_qp])*(-8*(1 + Utility::pow<2>(_alpha))*_nz*(2*_K1*(_nx + _ny) + 3*_K2*_nx*Utility::pow<2>(_nz))*std::cos(2*_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 128*_alpha*_K2*Utility::pow<2>(_ny)*(-3*Utility::pow<2>(_nx) + Utility::pow<2>(_ny))*Utility::pow<3>(std::sin(_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp])) + 
                (1 + Utility::pow<2>(_alpha))*_nz*(44*_K1*_nx + 4*_K1*_ny + 63*_K2*_nx*Utility::pow<2>(_ny) + 60*_K2*_nx*Utility::pow<2>(_nz) + 72*_K2*_nx*Utility::pow<2>(_ny)*std::cos(4*_azimuth_phi[_qp])*Utility::pow<4>(std::sin(_polar_theta[_qp])))) + 
             _nz*(-32*_alpha*_K2*_ny*std::cos(2*_polar_theta[_qp])*(-3*Utility::pow<2>(_ny) + 2*Utility::pow<2>(_nz) + 6*(-_nx + _ny)*(_nx + _ny)*std::cos(2*_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp]))) - 
                (1 + Utility::pow<2>(_alpha))*std::cos(3*_polar_theta[_qp])*(8*_K1*(3*_nx + _ny) + _K2*_nx*(261*Utility::pow<2>(_ny) - 68*Utility::pow<2>(_nz)) + 144*_K2*_nx*(-4*Utility::pow<2>(_ny) + Utility::pow<2>(_nz))*std::cos(2*_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp]))) + 
                _K2*(16*_alpha*_ny*(3*Utility::pow<2>(_nx) - 3*Utility::pow<2>(_ny) + 2*Utility::pow<2>(_nz))*std::cos(4*_polar_theta[_qp]) + 3*_nx*(5*(1 + Utility::pow<2>(_alpha))*(9*Utility::pow<2>(_ny) - 4*Utility::pow<2>(_nz))*std::cos(5*_polar_theta[_qp]) + 64*_alpha*_nx*_ny*std::cos(2*_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp])))))))/
        64. + Utility::pow<2>(std::cos(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp]))*(6*(1 + Utility::pow<2>(_alpha))*_K2*_ny*_nz*(-3*Utility::pow<2>(_nx) + Utility::pow<2>(_nz))*Utility::pow<3>(std::cos(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp]) + 
          3*(1 + Utility::pow<2>(_alpha))*_K2*_ny*_nz*(12*Utility::pow<2>(_nx) + Utility::pow<2>(_ny) - Utility::pow<2>(_ny)*std::cos(2*_azimuth_phi[_qp]))*std::cos(_polar_theta[_qp])*std::sin(_azimuth_phi[_qp])*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 
          6*_K2*_nz*Utility::pow<2>(std::cos(_polar_theta[_qp]))*(-(_alpha*_nx*(Utility::pow<2>(_nx) - 2*Utility::pow<2>(_ny))*std::sin(_azimuth_phi[_qp])) - (1 + Utility::pow<2>(_alpha))*_nz*(-3*Utility::pow<2>(_nx) + Utility::pow<2>(_nz))*std::sin(_polar_theta[_qp])) + 
          std::sin(_polar_theta[_qp])*((1 + Utility::pow<2>(_alpha))*_K1*(_nx*(_nx + _ny) - Utility::pow<2>(_nz)) + 6*_K2*_nz*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])*(_alpha*_nx*(Utility::pow<2>(_nx) - 2*Utility::pow<2>(_ny)) - (1 + Utility::pow<2>(_alpha))*(Utility::pow<2>(_nx) + Utility::pow<2>(_ny))*_nz*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]))) + 
          3*_K2*_ny*(-2*_alpha*_nx*(Utility::pow<2>(_nz) + (_nx - _ny)*(_nx + _ny)*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))) - (1 + Utility::pow<2>(_alpha))*_nz*std::sin(_azimuth_phi[_qp])*(2*Utility::pow<2>(_nz) - 3*Utility::pow<2>(_nx)*Utility::pow<2>(std::sin(_azimuth_phi[_qp])))*std::sin(_polar_theta[_qp]))*std::sin(2*_polar_theta[_qp])) + 
       (8*Utility::pow<2>(std::sin(_polar_theta[_qp]))*(_alpha*_K1*_nx*_nz*std::sin(_azimuth_phi[_qp]) - (1 + Utility::pow<2>(_alpha))*_K1*Utility::pow<2>(_nz)*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp]) + 6*_alpha*_K2*_nx*_ny*_nz*Utility::pow<2>(std::sin(_azimuth_phi[_qp]))*(2*_nx*std::cos(_azimuth_phi[_qp]) + _ny*std::sin(_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 
             (1 + Utility::pow<2>(_alpha))*_K2*Utility::pow<2>(_ny)*(2*(Utility::pow<2>(_ny) - 3*Utility::pow<2>(_nz))*Utility::pow<4>(std::sin(_azimuth_phi[_qp])) + 3*Utility::pow<2>(_nx)*Utility::pow<2>(std::sin(2*_azimuth_phi[_qp])))*Utility::pow<3>(std::sin(_polar_theta[_qp]))) + 
          8*(1 + Utility::pow<2>(_alpha))*_K1*_ny*_nz*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])*std::sin(2*_polar_theta[_qp]) + 3*_K2*(_alpha*_nx*Utility::pow<2>(_ny)*_nz*std::sin(3*_azimuth_phi[_qp]) - 
             (1 + Utility::pow<2>(_alpha))*(2*Utility::pow<2>(_nx)*Utility::pow<2>(_ny) + Utility::pow<4>(_ny) - (Utility::pow<2>(_nx) + 7*Utility::pow<2>(_ny))*Utility::pow<2>(_nz) + 2*Utility::pow<4>(_nz) - 8*Utility::pow<2>(_ny)*Utility::pow<2>(_nz)*Utility::pow<4>(std::sin(_azimuth_phi[_qp])))*std::sin(_polar_theta[_qp]))*Utility::pow<2>(std::sin(2*_polar_theta[_qp])) + 
          8*(1 + Utility::pow<2>(_alpha))*_K2*_ny*Utility::pow<3>(_nz)*(1.0/std::sin(_polar_theta[_qp]))*std::sin(_azimuth_phi[_qp])*Utility::pow<3>(std::sin(2*_polar_theta[_qp])))/8. + 
       _K2*Utility::pow<2>(_nx)*Utility::pow<4>(std::cos(_azimuth_phi[_qp]))*Utility::pow<3>(std::sin(_polar_theta[_qp]))*(4*_alpha*_nx*_ny*std::cos(_polar_theta[_qp]) + (1 + Utility::pow<2>(_alpha))*(-((Utility::pow<2>(_nx) - 3*Utility::pow<2>(_nz))*(1 + 3*std::cos(2*_polar_theta[_qp]))) + 9*_ny*_nz*std::sin(_azimuth_phi[_qp])*std::sin(2*_polar_theta[_qp])))))/(1 + Utility::pow<2>(_alpha)));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
