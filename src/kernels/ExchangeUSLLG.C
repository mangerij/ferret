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

#include "ExchangeUSLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", ExchangeUSLLG);

template<>
InputParameters validParams<ExchangeUSLLG>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for polar, 1 for azimuth)");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  params.addRequiredParam<Real>("alpha", "the damping coefficient in the LLG equation");
  params.addRequiredParam<Real>("g0", "g0");
  params.addRequiredParam<Real>("Ae", "Ae");
  params.addRequiredParam<Real>("Ms", "Ms");
  return params;
}

ExchangeUSLLG::ExchangeUSLLG(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _azimuth_phi_var(coupled("azimuth_phi")),
  _polar_theta_var(coupled("polar_theta")),
  _azimuth_phi(coupledValue("azimuth_phi")),
  _polar_theta(coupledValue("polar_theta")),
  _azimuth_phi_grad(coupledGradient("azimuth_phi")),
  _polar_theta_grad(coupledGradient("polar_theta")),
  _alpha(getParam<Real>("alpha")),
  _g0(getParam<Real>("g0")),
  _Ae(getParam<Real>("Ae")),
  _Ms(getParam<Real>("Ms"))
{
}

Real
ExchangeUSLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return (-2*_Ae*_g0*(_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2) + _alpha*(Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2)))*_test[_i][_qp]*std::cos(_polar_theta[_qp]) - (_azimuth_phi_grad[_qp](0)*_polar_theta_grad[_qp](0) + _azimuth_phi_grad[_qp](1)*_polar_theta_grad[_qp](1) + _azimuth_phi_grad[_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp]*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])) + 
       _alpha*(_grad_test[_i][_qp](0)*_polar_theta_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_theta_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_theta_grad[_qp](2))*(1.0/std::sin(_polar_theta[_qp])))*std::sqrt(Utility::pow<2>(std::sin(_polar_theta[_qp]))))/((1 + Utility::pow<2>(_alpha))*_Ms);
  }
  else if (_component == 1)
  {
    return -(_Ae*_g0*(1.0/std::sin(_polar_theta[_qp]))*(4*(_grad_test[_i][_qp](0)*_polar_theta_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_theta_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_theta_grad[_qp](2))*std::cos(2*_azimuth_phi[_qp])*Utility::pow<2>(std::cos(_polar_theta[_qp])) + 
        (Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_test[_i][_qp]*(2 - 6*std::cos(2*_azimuth_phi[_qp]) - 2*std::cos(2*_polar_theta[_qp]) + std::cos(2*(-_azimuth_phi[_qp] + _polar_theta[_qp])) + std::cos(2*(_azimuth_phi[_qp] + _polar_theta[_qp])))*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])) + 
        4*std::cos(_polar_theta[_qp])*(-2*_alpha*(_azimuth_phi_grad[_qp](0)*_polar_theta_grad[_qp](0) + _azimuth_phi_grad[_qp](1)*_polar_theta_grad[_qp](1) + _azimuth_phi_grad[_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp] + ((Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_test[_i][_qp] + (_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*std::sin(2*_azimuth_phi[_qp]))*
            std::sin(_polar_theta[_qp])) + 2*(2*_alpha*(_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*std::sin(_polar_theta[_qp]) + 
           2*(_grad_test[_i][_qp](0)*_polar_theta_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_theta_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_theta_grad[_qp](2) - (_azimuth_phi_grad[_qp](0)*_polar_theta_grad[_qp](0) + _azimuth_phi_grad[_qp](1)*_polar_theta_grad[_qp](1) + _azimuth_phi_grad[_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp]*std::sin(2*_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 
           (3*Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + 3*Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + 3*Utility::pow<2>(_azimuth_phi_grad[_qp](2)) + Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_test[_i][_qp]*std::cos(2*_azimuth_phi[_qp])*std::sin(2*_polar_theta[_qp]))))/(2.*(1 + Utility::pow<2>(_alpha))*_Ms);
  }
  else
    return 0.0;
}

Real
ExchangeUSLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return (-2*_Ae*_g0*std::sin(_polar_theta[_qp])*(_alpha*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)) + ((_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*_phi[_j][_qp] - (_azimuth_phi_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_phi[_j][_qp](2))*_test[_i][_qp])*std::cos(_polar_theta[_qp]) + 
       _alpha*(Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + Utility::pow<2>(_azimuth_phi_grad[_qp](2)))*_phi[_j][_qp]*_test[_i][_qp]*std::cos(2*_polar_theta[_qp]) + (_azimuth_phi_grad[_qp](0)*_polar_theta_grad[_qp](0) + _azimuth_phi_grad[_qp](1)*_polar_theta_grad[_qp](1) + _azimuth_phi_grad[_qp](2)*_polar_theta_grad[_qp](2))*_phi[_j][_qp]*_test[_i][_qp]*std::sin(_polar_theta[_qp])))/((1 + Utility::pow<2>(_alpha))*_Ms*std::sqrt(Utility::pow<2>(std::sin(_polar_theta[_qp]))));
  }
  else if (_component == 1)
  {
    return (-2*_Ae*_g0*(1.0/std::sin(_polar_theta[_qp]))*(-2*(_grad_test[_i][_qp](0)*_polar_theta_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_theta_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_theta_grad[_qp](2))*_phi[_j][_qp]*Utility::pow<2>(std::cos(_polar_theta[_qp]))*std::sin(2*_azimuth_phi[_qp]) - 
       (Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_phi[_j][_qp]*_test[_i][_qp]*(-3 + std::cos(2*_polar_theta[_qp]))*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*std::sin(2*_azimuth_phi[_qp]) + _alpha*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*std::sin(_polar_theta[_qp]) + 
       (-2*(_azimuth_phi_grad[_qp](0)*_polar_theta_grad[_qp](0) + _azimuth_phi_grad[_qp](1)*_polar_theta_grad[_qp](1) + _azimuth_phi_grad[_qp](2)*_polar_theta_grad[_qp](2))*_phi[_j][_qp]*_test[_i][_qp]*std::cos(2*_azimuth_phi[_qp]) - (_grad_phi[_j][_qp](0)*_polar_theta_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_theta_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp]*std::sin(2*_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 
       std::cos(_polar_theta[_qp])*(-2*_alpha*(_grad_phi[_j][_qp](0)*_polar_theta_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_theta_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp] + 2*(_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*_phi[_j][_qp]*std::cos(2*_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]) + 
          (_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*std::sin(2*_azimuth_phi[_qp])*std::sin(_polar_theta[_qp])) + 3*(_azimuth_phi_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_phi[_j][_qp](2))*_test[_i][_qp]*std::cos(2*_azimuth_phi[_qp])*std::sin(2*_polar_theta[_qp]) - 
       (3*Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + 3*Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + 3*Utility::pow<2>(_azimuth_phi_grad[_qp](2)) + Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_phi[_j][_qp]*_test[_i][_qp]*std::sin(2*_azimuth_phi[_qp])*std::sin(2*_polar_theta[_qp])))/((1 + Utility::pow<2>(_alpha))*_Ms);
  }
  else
    return 0.0;
}

Real
ExchangeUSLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuth_phi_var)
    {
      return (-2*_Ae*_g0*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2) + 2*_alpha*(_azimuth_phi_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_phi[_j][_qp](2))*_test[_i][_qp]*std::cos(_polar_theta[_qp]) - (_grad_phi[_j][_qp](0)*_polar_theta_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_theta_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp]*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])))*
     std::sqrt(Utility::pow<2>(std::sin(_polar_theta[_qp]))))/((1 + Utility::pow<2>(_alpha))*_Ms);
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
      return (_Ae*_g0*(1.0/std::sin(_polar_theta[_qp]))*(-4*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*std::cos(2*_azimuth_phi[_qp])*Utility::pow<2>(std::cos(_polar_theta[_qp])) - 
       2*(_grad_phi[_j][_qp](0)*_polar_theta_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_theta_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp]*(2 - 6*std::cos(2*_azimuth_phi[_qp]) - 2*std::cos(2*_polar_theta[_qp]) + std::cos(2*(-_azimuth_phi[_qp] + _polar_theta[_qp])) + std::cos(2*(_azimuth_phi[_qp] + _polar_theta[_qp])))*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])) + 
       (Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_phi[_j][_qp]*_test[_i][_qp]*(2 - 6*std::cos(2*_azimuth_phi[_qp]) - 2*std::cos(2*_polar_theta[_qp]) + std::cos(2*(-_azimuth_phi[_qp] + _polar_theta[_qp])) + std::cos(2*(_azimuth_phi[_qp] + _polar_theta[_qp])))*Utility::pow<2>((1.0/std::sin(_polar_theta[_qp]))) - 
       16*(Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>(std::cos(_polar_theta[_qp]))*Utility::pow<2>(std::sin(_azimuth_phi[_qp])) - 
       4*std::cos(_polar_theta[_qp])*(-2*_alpha*(_azimuth_phi_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_phi[_j][_qp](2))*_test[_i][_qp] + _phi[_j][_qp]*std::cos(_polar_theta[_qp])*
           ((Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_test[_i][_qp] + (_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*std::sin(2*_azimuth_phi[_qp])) + 2*(_grad_phi[_j][_qp](0)*_polar_theta_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_theta_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp]*std::sin(_polar_theta[_qp])) + 
       4*_phi[_j][_qp]*std::sin(_polar_theta[_qp])*(-2*_alpha*(_azimuth_phi_grad[_qp](0)*_polar_theta_grad[_qp](0) + _azimuth_phi_grad[_qp](1)*_polar_theta_grad[_qp](1) + _azimuth_phi_grad[_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp] + 
          ((Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_test[_i][_qp] + (_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*std::sin(2*_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp])) + 
       4*(_grad_test[_i][_qp](0)*_polar_theta_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_theta_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_theta_grad[_qp](2))*_phi[_j][_qp]*std::cos(2*_azimuth_phi[_qp])*std::sin(2*_polar_theta[_qp]) - 
       4*(_alpha*(_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*_phi[_j][_qp]*std::cos(_polar_theta[_qp]) + (3*Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + 3*Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + 3*Utility::pow<2>(_azimuth_phi_grad[_qp](2)) + Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_phi[_j][_qp]*_test[_i][_qp]*
           std::cos(2*_azimuth_phi[_qp])*std::cos(2*_polar_theta[_qp]) + (_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2) - (_azimuth_phi_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_phi[_j][_qp](2))*_test[_i][_qp]*std::sin(2*_azimuth_phi[_qp]))*Utility::pow<2>(std::sin(_polar_theta[_qp])) + 
          (_grad_phi[_j][_qp](0)*_polar_theta_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_theta_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp]*std::cos(2*_azimuth_phi[_qp])*std::sin(2*_polar_theta[_qp]) + 
          _phi[_j][_qp]*(_grad_test[_i][_qp](0)*_polar_theta_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_theta_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_theta_grad[_qp](2) - (_azimuth_phi_grad[_qp](0)*_polar_theta_grad[_qp](0) + _azimuth_phi_grad[_qp](1)*_polar_theta_grad[_qp](1) + _azimuth_phi_grad[_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp]*std::sin(2*_azimuth_phi[_qp]))*std::sin(2*_polar_theta[_qp])) + 
       _phi[_j][_qp]*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp]))*(4*(_grad_test[_i][_qp](0)*_polar_theta_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_theta_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_theta_grad[_qp](2))*std::cos(2*_azimuth_phi[_qp])*Utility::pow<2>(std::cos(_polar_theta[_qp])) + 
          (Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_test[_i][_qp]*(2 - 6*std::cos(2*_azimuth_phi[_qp]) - 2*std::cos(2*_polar_theta[_qp]) + std::cos(2*(-_azimuth_phi[_qp] + _polar_theta[_qp])) + std::cos(2*(_azimuth_phi[_qp] + _polar_theta[_qp])))*(std::cos(_polar_theta[_qp])/std::sin(_polar_theta[_qp])) + 
          4*std::cos(_polar_theta[_qp])*(-2*_alpha*(_azimuth_phi_grad[_qp](0)*_polar_theta_grad[_qp](0) + _azimuth_phi_grad[_qp](1)*_polar_theta_grad[_qp](1) + _azimuth_phi_grad[_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp] + 
             ((Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_test[_i][_qp] + (_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*std::sin(2*_azimuth_phi[_qp]))*std::sin(_polar_theta[_qp])) + 
          2*(2*_alpha*(_azimuth_phi_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuth_phi_grad[_qp](1)*_grad_test[_i][_qp](1) + _azimuth_phi_grad[_qp](2)*_grad_test[_i][_qp](2))*std::sin(_polar_theta[_qp]) + 2*(_grad_test[_i][_qp](0)*_polar_theta_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_theta_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_theta_grad[_qp](2) - (_azimuth_phi_grad[_qp](0)*_polar_theta_grad[_qp](0) + _azimuth_phi_grad[_qp](1)*_polar_theta_grad[_qp](1) + _azimuth_phi_grad[_qp](2)*_polar_theta_grad[_qp](2))*_test[_i][_qp]*std::sin(2*_azimuth_phi[_qp]))*
              Utility::pow<2>(std::sin(_polar_theta[_qp])) + (3*Utility::pow<2>(_azimuth_phi_grad[_qp](0)) + 3*Utility::pow<2>(_azimuth_phi_grad[_qp](1)) + 3*Utility::pow<2>(_azimuth_phi_grad[_qp](2)) + Utility::pow<2>(_polar_theta_grad[_qp](0)) + Utility::pow<2>(_polar_theta_grad[_qp](1)) + Utility::pow<2>(_polar_theta_grad[_qp](2)))*_test[_i][_qp]*std::cos(2*_azimuth_phi[_qp])*std::sin(2*_polar_theta[_qp])))))/(2.*(1 + Utility::pow<2>(_alpha))*_Ms);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
