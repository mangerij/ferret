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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ExchangeUSLL.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", ExchangeUSLL);

template<>
InputParameters validParams<ExchangeUSLL>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for polar, 1 for azimuth)");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  return params;
}

ExchangeUSLL::ExchangeUSLL(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _azimuth_phi_var(coupled("azimuth_phi")),
  _polar_theta_var(coupled("polar_theta")),
  _azimuth_phi(coupledValue("azimuth_phi")),
  _polar_theta(coupledValue("polar_theta")),
  _azimuth_phi_grad(coupledGradient("azimuth_phi")),
  _polar_theta_grad(coupledGradient("polar_theta")),
  _alpha(getMaterialProperty<Real>("alpha")),
  _g0(getMaterialProperty<Real>("g0")),
  _Ae(getMaterialProperty<Real>("Ae")),
  _Ms(getMaterialProperty<Real>("Ms"))
{
}

Real
ExchangeUSLL::computeQpResidual()
{
  if (_component == 0)
  {
    return (-2*_Ae[_qp]*_g0[_qp]*(std::sin(_polar_theta[_qp])*(_alpha[_qp]*(_grad_test[_i][_qp](0)*_polar_theta_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_theta_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_theta_grad[_qp](2)) + (_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*std::sin(_polar_theta[_qp])) + 
       _alpha[_qp]*_test[_i][_qp]*std::cos(_polar_theta[_qp])*(Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)) + (Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2)))*Utility::pow<2>(std::sin(_polar_theta[_qp])))))/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 1)
  {
    return (-2*_Ae[_qp]*_g0[_qp]*(_grad_test[_i][_qp](0)*_polar_theta_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_theta_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_theta_grad[_qp](2) - _alpha[_qp]*(_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*std::sin(_polar_theta[_qp]) + 
       _test[_i][_qp]*std::cos(_polar_theta[_qp])*(_alpha[_qp]*(_azimuth_phi_grad[_qp](0)*_polar_theta_grad[_qp](0) + _azimuth_phi_grad[_qp](1)*_polar_theta_grad[_qp](1) + _azimuth_phi_grad[_qp](2)*_polar_theta_grad[_qp](2)) + (Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2)))*std::sin(_polar_theta[_qp]))))/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else
    return 0.0;
}

Real
ExchangeUSLL::computeQpJacobian()
{
  if (_component == 0)
  {
    return (_Ae[_qp]*_g0[_qp]*(-4*_alpha[_qp]*((_grad_test[_i][_qp](0)*_polar_theta_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_theta_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_theta_grad[_qp](2))*_phi[_j][_qp] + 2*(_grad_phi[_j][_qp](0)*_polar_theta_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_theta_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp])*std::cos(_polar_theta[_qp]) + 
       _alpha[_qp]*(-4*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)) + (Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2)) + 4*(Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2))))*_phi[_j][_qp]*_test[_i][_qp])*std::sin(_polar_theta[_qp]) - 
       4*(_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*_phi[_j][_qp]*std::sin(2*_polar_theta[_qp]) - 3*_alpha[_qp]*(Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2)))*_phi[_j][_qp]*_test[_i][_qp]*std::sin(3*_polar_theta[_qp])))/(2.*(1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 1)
  {
    return (-2*_Ae[_qp]*_g0[_qp]*(-(_alpha[_qp]*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*std::sin(_polar_theta[_qp])) + _test[_i][_qp]*std::cos(_polar_theta[_qp])*
        (_alpha[_qp]*(_grad_phi[_j][_qp](0)*_polar_theta_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_theta_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_theta_grad[_qp](2)) + 2*(_azimuth_phi_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_phi[_j][_qp](2))*std::sin(_polar_theta[_qp]))))/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else
    return 0.0;
}

Real
ExchangeUSLL::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuth_phi_var)
    {
      return (-2*_Ae[_qp]*_g0[_qp]*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2) + 2*_alpha[_qp]*(_azimuth_phi_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_phi[_j][_qp](2))*_test[_i][_qp]*std::cos(_polar_theta[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])))/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
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
      return (2*_Ae[_qp]*_g0[_qp]*(-(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)) - _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) - _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2) + _alpha[_qp]*((_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*_phi[_j][_qp] - (_azimuth_phi_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_phi[_j][_qp](2))*_test[_i][_qp])*std::cos(_polar_theta[_qp]) - 
       (Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2)))*_phi[_j][_qp]*_test[_i][_qp]*std::cos(2*_polar_theta[_qp]) + _alpha[_qp]*(_azimuth_phi_grad[_qp](0)*_polar_theta_grad[_qp](0) + _azimuth_phi_grad[_qp](1)*_polar_theta_grad[_qp](1) + _azimuth_phi_grad[_qp](2)*_polar_theta_grad[_qp](2))*_phi[_j][_qp]*_test[_i][_qp]*std::sin(_polar_theta[_qp])))/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
