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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "MasterExchangeUSLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MasterExchangeUSLLG);

InputParameters MasterExchangeUSLLG::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to the magnetic exchange energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_th", "The polar angle of the constrained magnetic vector");
  params.addRequiredCoupledVar("azimuthal_ph", "The azimuthal angle of the constrained magnetic vector");
  return params;
}

MasterExchangeUSLLG::MasterExchangeUSLLG(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _polar_th_var(coupled("polar_th")),
  _azimuthal_ph_var(coupled("azimuthal_ph")),
  _polar_th(coupledValue("polar_th")),
  _azimuthal_ph(coupledValue("azimuthal_ph")),
  _polar_th_grad(coupledGradient("polar_th")),
  _azimuthal_ph_grad(coupledGradient("azimuthal_ph")),
  _alpha(getMaterialProperty<Real>("alpha")),
  _g0(getMaterialProperty<Real>("g0")),
  _Ae(getMaterialProperty<Real>("Ae")),
  _Ms(getMaterialProperty<Real>("Ms"))
{
}

Real
MasterExchangeUSLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return (_Ae[_qp]*_g0[_qp]*(-2.0*_alpha[_qp]*(_grad_test[_i][_qp](0)*_polar_th_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_th_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_th_grad[_qp](2)) + 
       (4.0*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2)) - 2.0*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](2))*_test[_i][_qp])*std::cos(_polar_th[_qp]) - 
       2.0*_alpha[_qp]*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)))*_polar_th_grad[_qp](2)*_test[_i][_qp]*std::cos(2.0*_polar_th[_qp]) - 
       2.0*(_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2) + 2.0*_polar_th_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp])*std::sin(_polar_th[_qp]) - 
       _alpha[_qp]*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)))*_grad_test[_i][_qp](2)*std::sin(2.0*_polar_th[_qp])))/((1.0 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 1)
  {
    return (2.0*_Ae[_qp]*_g0[_qp]*(_grad_test[_i][_qp](0)*_polar_th_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_th_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_th_grad[_qp](2) + _polar_th_grad[_qp](2)*(Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
       Utility::pow<2>(std::cos(_azimuthal_ph[_qp]))*(_polar_th_grad[_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] - 
          (_grad_test[_i][_qp](0)*_polar_th_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_th_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_th_grad[_qp](2) + 4.0*_azimuthal_ph_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp])*Utility::pow<2>((std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))) - 
          _alpha[_qp]*(_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2))*(1.0/std::sin(_polar_th[_qp])) - 
          2.0*_alpha[_qp]*_polar_th_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*Utility::pow<3>((1.0/std::sin(_polar_th[_qp]))) + 
          (std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*(-(_grad_test[_i][_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))) + 
             (2.0*Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + 4.0*_azimuthal_ph_grad[_qp](1)*_azimuthal_ph_grad[_qp](2) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
             2.0*_alpha[_qp]*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*(1.0/std::sin(_polar_th[_qp])) + 
             (Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp]*Utility::pow<2>((1.0/std::sin(_polar_th[_qp]))))) + 
       2.0*std::cos(_azimuthal_ph[_qp])*(2.0*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2)) - 
          (3.0*_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*(2.0*_polar_th_grad[_qp](1) + _polar_th_grad[_qp](2)) + _azimuthal_ph_grad[_qp](2)*(_polar_th_grad[_qp](1) + 2.0*_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
          2.0*(Utility::pow<3>(_azimuthal_ph_grad[_qp](2)) + (_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1))*_polar_th_grad[_qp](2) + 
             _azimuthal_ph_grad[_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + 2.0*Utility::pow<2>(_polar_th_grad[_qp](2))))*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])) + 
          2.0*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*(-_grad_test[_i][_qp](2) + _test[_i][_qp] + _polar_th_grad[_qp](2)*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*Utility::pow<2>((1.0/std::sin(_polar_th[_qp]))))*std::sin(_azimuthal_ph[_qp]) + 
       (_grad_test[_i][_qp](0)*_polar_th_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_th_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_th_grad[_qp](2) + 4.0*_azimuthal_ph_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp])*Utility::pow<2>((std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*
        Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) - (_polar_th_grad[_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*
           _test[_i][_qp] + _alpha[_qp]*(1.0/std::sin(_polar_th[_qp]))*(_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2) + 
             2.0*_polar_th_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*Utility::pow<2>((1.0/std::sin(_polar_th[_qp])))))*Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) + 
       (std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*(-((Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*(_grad_test[_i][_qp](2) - _test[_i][_qp])) + 
          (_grad_test[_i][_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2))) - 
             (2.0*Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + 4.0*_azimuthal_ph_grad[_qp](1)*_azimuthal_ph_grad[_qp](2) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
             2.0*_alpha[_qp]*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*(1.0/std::sin(_polar_th[_qp])) - 
             (Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp]*Utility::pow<2>((1.0/std::sin(_polar_th[_qp]))))*Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) + 
          (_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2))*std::sin(2.0*_azimuthal_ph[_qp])))*std::sin(_polar_th[_qp]))/((1.0 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else
    return 0.0;
}

Real
MasterExchangeUSLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return (2.0*_Ae[_qp]*_g0[_qp]*(-(_alpha[_qp]*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))) + 
       (2.0*_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](2) + 2.0*_azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](2) + 2.0*_azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2) - _azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0)*_phi[_j][_qp] - _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1)*_phi[_j][_qp] - 
          _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2)*_phi[_j][_qp] - (_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](2) + 2.0*_polar_th_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_phi[_j][_qp])*_test[_i][_qp])*
        std::cos(_polar_th[_qp]) - _alpha[_qp]*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)))*(_grad_test[_i][_qp](2)*_phi[_j][_qp] + _grad_phi[_j][_qp](2)*_test[_i][_qp])*std::cos(2.0*_polar_th[_qp]) - 
       (2.0*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_phi[_j][_qp] + 
          2.0*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](2) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1)*_polar_th_grad[_qp](2) + 2.0*_azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp] - 
          (_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](2))*_phi[_j][_qp]*_test[_i][_qp])*std::sin(_polar_th[_qp]) + 
       2.0*_alpha[_qp]*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)))*_polar_th_grad[_qp](2)*_phi[_j][_qp]*_test[_i][_qp]*std::sin(2.0*_polar_th[_qp])))/((1.0 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 1)
  {
    return (2.0*_Ae[_qp]*_g0[_qp]*(2.0*_phi[_j][_qp]*Utility::pow<2>(std::cos(_azimuthal_ph[_qp]))*(2.0*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2)) - 
          (3.0*_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*(2.0*_polar_th_grad[_qp](1) + _polar_th_grad[_qp](2)) + _azimuthal_ph_grad[_qp](2)*(_polar_th_grad[_qp](1) + 2.0*_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
          2.0*(Utility::pow<3>(_azimuthal_ph_grad[_qp](2)) + (_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1))*_polar_th_grad[_qp](2) + 
             _azimuthal_ph_grad[_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + 2.0*Utility::pow<2>(_polar_th_grad[_qp](2))))*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])) + 
          2.0*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*(-_grad_test[_i][_qp](2) + _test[_i][_qp] + _polar_th_grad[_qp](2)*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*Utility::pow<2>((1.0/std::sin(_polar_th[_qp])))) + 
       Utility::pow<2>(std::cos(_azimuthal_ph[_qp]))*(2.0*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2))*_polar_th_grad[_qp](2)*_test[_i][_qp] + 
          2.0*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*(-((_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2))*_grad_test[_i][_qp](2)) + 2.0*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](2))*_test[_i][_qp] - 
             2.0*(_azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](1) + 2.0*_azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))) - 
          _alpha[_qp]*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](2) - 2.0*_grad_test[_i][_qp](2)*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*
           (1.0/std::sin(_polar_th[_qp])) - 2.0*_alpha[_qp]*_polar_th_grad[_qp](2)*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*Utility::pow<3>((1.0/std::sin(_polar_th[_qp])))) + 
       2.0*_phi[_j][_qp]*(_grad_test[_i][_qp](0)*_polar_th_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_th_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_th_grad[_qp](2) + 4.0*_azimuthal_ph_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp])*std::cos(_azimuthal_ph[_qp])*
        Utility::pow<2>((std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*std::sin(_azimuthal_ph[_qp]) + 2.0*std::cos(_azimuthal_ph[_qp])*(2.0*_grad_test[_i][_qp](2)*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2)) - 
          (3.0*_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*(2.0*_polar_th_grad[_qp](1) + _polar_th_grad[_qp](2)) + _grad_phi[_j][_qp](2)*(_polar_th_grad[_qp](1) + 2.0*_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
          2.0*(2.0*_azimuthal_ph_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1)) + 3.0*Utility::pow<2>(_azimuthal_ph_grad[_qp](2))*_grad_phi[_j][_qp](2) + (_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1))*_polar_th_grad[_qp](2) + 
             _grad_phi[_j][_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + 2.0*Utility::pow<2>(_polar_th_grad[_qp](2))))*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])) + 
          2.0*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*(-_grad_test[_i][_qp](2) + _test[_i][_qp] + _polar_th_grad[_qp](2)*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*Utility::pow<2>((1.0/std::sin(_polar_th[_qp]))))*std::sin(_azimuthal_ph[_qp]) - 
       2.0*_phi[_j][_qp]*std::cos(_azimuthal_ph[_qp])*(_polar_th_grad[_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
          _alpha[_qp]*(1.0/std::sin(_polar_th[_qp]))*(_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2) + 
             2.0*_polar_th_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*Utility::pow<2>((1.0/std::sin(_polar_th[_qp])))))*std::sin(_azimuthal_ph[_qp]) - 
       2.0*_phi[_j][_qp]*std::cos(_azimuthal_ph[_qp])*(_polar_th_grad[_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] - 
          (_grad_test[_i][_qp](0)*_polar_th_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_th_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_th_grad[_qp](2) + 4.0*_azimuthal_ph_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp])*Utility::pow<2>((std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))) - 
          _alpha[_qp]*(_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2))*(1.0/std::sin(_polar_th[_qp])) - 
          2.0*_alpha[_qp]*_polar_th_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*Utility::pow<3>((1.0/std::sin(_polar_th[_qp]))) + 
          (std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*(-(_grad_test[_i][_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))) + 
             (2.0*Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + 4.0*_azimuthal_ph_grad[_qp](1)*_azimuthal_ph_grad[_qp](2) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
             2.0*_alpha[_qp]*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*(1.0/std::sin(_polar_th[_qp])) + 
             (Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp]*Utility::pow<2>((1.0/std::sin(_polar_th[_qp])))))*std::sin(_azimuthal_ph[_qp]) + 
       4.0*(_azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](1) + 2.0*_azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*Utility::pow<2>((std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*
        Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) - 2.0*_phi[_j][_qp]*(2.0*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2)) - 
          (3.0*_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*(2.0*_polar_th_grad[_qp](1) + _polar_th_grad[_qp](2)) + _azimuthal_ph_grad[_qp](2)*(_polar_th_grad[_qp](1) + 2.0*_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
          2.0*(Utility::pow<3>(_azimuthal_ph_grad[_qp](2)) + (_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1))*_polar_th_grad[_qp](2) + 
             _azimuthal_ph_grad[_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + 2.0*Utility::pow<2>(_polar_th_grad[_qp](2))))*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])) + 
          2.0*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*(-_grad_test[_i][_qp](2) + _test[_i][_qp] + _polar_th_grad[_qp](2)*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*Utility::pow<2>((1.0/std::sin(_polar_th[_qp]))))*Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) - 
       (2.0*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2))*_polar_th_grad[_qp](2)*_test[_i][_qp] + 
          _alpha[_qp]*(1.0/std::sin(_polar_th[_qp]))*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](2) + 
             2.0*_polar_th_grad[_qp](2)*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*Utility::pow<2>((1.0/std::sin(_polar_th[_qp])))))*Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) + 
       (std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*(2.0*(_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2))*_phi[_j][_qp]*std::cos(2.0*_azimuthal_ph[_qp]) + 
          2.0*((_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2))*_grad_test[_i][_qp](2) - 2.0*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](2))*_test[_i][_qp] + 
             _alpha[_qp]*_grad_test[_i][_qp](2)*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*(1.0/std::sin(_polar_th[_qp])))*Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) + 
          (_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](1) + _grad_test[_i][_qp](2)*(_grad_phi[_j][_qp](1) + 
                (Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_phi[_j][_qp]) - 
             (2.0*Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + 4.0*_azimuthal_ph_grad[_qp](1)*_azimuthal_ph_grad[_qp](2) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_phi[_j][_qp]*_test[_i][_qp] + 
             _phi[_j][_qp]*(1.0/std::sin(_polar_th[_qp]))*(2.0*_alpha[_qp]*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2)) - 
                (Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp]*(1.0/std::sin(_polar_th[_qp]))))*std::sin(2.0*_azimuthal_ph[_qp])))*std::sin(_polar_th[_qp]))/((1.0 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else
    return 0.0;
}

Real
MasterExchangeUSLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _azimuthal_ph_var)
    {
      return (2.0*_Ae[_qp]*_g0[_qp]*((2.0*_grad_test[_i][_qp](2)*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2)) - (_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](2))*_test[_i][_qp])*std::cos(_polar_th[_qp]) - 
       2.0*_alpha[_qp]*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2))*_polar_th_grad[_qp](2)*_test[_i][_qp]*std::cos(2.0*_polar_th[_qp]) - 
       (_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](2) + 2.0*_polar_th_grad[_qp](2)*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp])*std::sin(_polar_th[_qp]) - 
       _alpha[_qp]*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2))*_grad_test[_i][_qp](2)*std::sin(2.0*_polar_th[_qp])))/((1.0 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _polar_th_var)
    {
      return (2.0*_Ae[_qp]*_g0[_qp]*(_phi[_j][_qp]*std::cos(_polar_th[_qp])*(_grad_test[_i][_qp](0)*_polar_th_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_th_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_th_grad[_qp](2) + _polar_th_grad[_qp](2)*(Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
          Utility::pow<2>(std::cos(_azimuthal_ph[_qp]))*(_polar_th_grad[_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*
              _test[_i][_qp] - (_grad_test[_i][_qp](0)*_polar_th_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_th_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_th_grad[_qp](2) + 4.0*_azimuthal_ph_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp])*
              Utility::pow<2>((std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))) - _alpha[_qp]*(_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2))*(1.0/std::sin(_polar_th[_qp])) - 
             2.0*_alpha[_qp]*_polar_th_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*Utility::pow<3>((1.0/std::sin(_polar_th[_qp]))) + 
             (std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*(-(_grad_test[_i][_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))) + 
                (2.0*Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + 4.0*_azimuthal_ph_grad[_qp](1)*_azimuthal_ph_grad[_qp](2) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
                (1.0/std::sin(_polar_th[_qp]))*(2.0*_alpha[_qp]*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2)) + 
                   (Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp]*(1.0/std::sin(_polar_th[_qp]))))) + 
          2.0*std::cos(_azimuthal_ph[_qp])*(2.0*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2)) - 
             (3.0*_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*(2.0*_polar_th_grad[_qp](1) + _polar_th_grad[_qp](2)) + _azimuthal_ph_grad[_qp](2)*(_polar_th_grad[_qp](1) + 2.0*_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
             2.0*(Utility::pow<3>(_azimuthal_ph_grad[_qp](2)) + (_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1))*_polar_th_grad[_qp](2) + 
                _azimuthal_ph_grad[_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + 2.0*Utility::pow<2>(_polar_th_grad[_qp](2))))*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])) + 
             2.0*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*(-_grad_test[_i][_qp](2) + _test[_i][_qp] + _polar_th_grad[_qp](2)*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*Utility::pow<2>((1.0/std::sin(_polar_th[_qp]))))*std::sin(_azimuthal_ph[_qp]) + 
          (_grad_test[_i][_qp](0)*_polar_th_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_th_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_th_grad[_qp](2) + 4.0*_azimuthal_ph_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp])*Utility::pow<2>((std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*
           Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) - (_polar_th_grad[_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*
              _test[_i][_qp] + _alpha[_qp]*(1.0/std::sin(_polar_th[_qp]))*(_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2) + 
                2.0*_polar_th_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*Utility::pow<2>((1.0/std::sin(_polar_th[_qp])))))*Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) + 
          (std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*(-((Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*(_grad_test[_i][_qp](2) - _test[_i][_qp])) + 
             (_grad_test[_i][_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2))) - 
                (2.0*Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + 4.0*_azimuthal_ph_grad[_qp](1)*_azimuthal_ph_grad[_qp](2) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
                (1.0/std::sin(_polar_th[_qp]))*(2.0*_alpha[_qp]*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2)) - 
                   (Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp]*(1.0/std::sin(_polar_th[_qp]))))*Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) + 
             (_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2))*std::sin(2.0*_azimuthal_ph[_qp]))) + 
       std::sin(_polar_th[_qp])*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2) + 2.0*_polar_th_grad[_qp](2)*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp] + 
          _grad_phi[_j][_qp](2)*(Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
          Utility::pow<2>(std::cos(_azimuthal_ph[_qp]))*(2.0*_polar_th_grad[_qp](2)*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp] + 
             _grad_phi[_j][_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] - 
             (_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2) + 4.0*_azimuthal_ph_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2))*_test[_i][_qp])*
              Utility::pow<2>((std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))) + _alpha[_qp]*(_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2))*_phi[_j][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*(1.0/std::sin(_polar_th[_qp])) + 
             2.0*_phi[_j][_qp]*(_grad_test[_i][_qp](0)*_polar_th_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_th_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_th_grad[_qp](2) + 4.0*_azimuthal_ph_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp])*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*
              Utility::pow<2>((1.0/std::sin(_polar_th[_qp]))) - 2.0*_alpha[_qp]*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2))*_polar_th_grad[_qp](2)*_test[_i][_qp]*Utility::pow<3>((1.0/std::sin(_polar_th[_qp]))) - 
             2.0*_alpha[_qp]*_grad_phi[_j][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*Utility::pow<3>((1.0/std::sin(_polar_th[_qp]))) + 
             6*_alpha[_qp]*_polar_th_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_phi[_j][_qp]*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*Utility::pow<3>((1.0/std::sin(_polar_th[_qp]))) + 
             2.0*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*(-((_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*(_grad_test[_i][_qp](2) - _test[_i][_qp])) + 
                _alpha[_qp]*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2) - (_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_phi[_j][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*(1.0/std::sin(_polar_th[_qp])) + 
                _test[_i][_qp]*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2) - (Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_phi[_j][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*
                 Utility::pow<2>((1.0/std::sin(_polar_th[_qp])))) - _phi[_j][_qp]*Utility::pow<2>((1.0/std::sin(_polar_th[_qp])))*
              (-(_grad_test[_i][_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))) + 
                (2.0*Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + 4.0*_azimuthal_ph_grad[_qp](1)*_azimuthal_ph_grad[_qp](2) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
                (1.0/std::sin(_polar_th[_qp]))*(2.0*_alpha[_qp]*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2)) + 
                   (Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp]*(1.0/std::sin(_polar_th[_qp]))))) + 
          (_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2) + 4.0*_azimuthal_ph_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2))*_test[_i][_qp])*Utility::pow<2>((std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*
           Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) - 2.0*_phi[_j][_qp]*(_grad_test[_i][_qp](0)*_polar_th_grad[_qp](0) + _grad_test[_i][_qp](1)*_polar_th_grad[_qp](1) + _grad_test[_i][_qp](2)*_polar_th_grad[_qp](2) + 4.0*_azimuthal_ph_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp])*
           (std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*Utility::pow<2>((1.0/std::sin(_polar_th[_qp])))*Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) - (-(_alpha[_qp]*_phi[_j][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*(1.0/std::sin(_polar_th[_qp]))*
                (_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2) + 6*_polar_th_grad[_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*Utility::pow<2>((1.0/std::sin(_polar_th[_qp]))))) + 
             _test[_i][_qp]*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0))*_grad_phi[_j][_qp](2) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1))*_grad_phi[_j][_qp](2) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2))*_grad_phi[_j][_qp](2) + _grad_phi[_j][_qp](2)*Utility::pow<2>(_polar_th_grad[_qp](0)) + _grad_phi[_j][_qp](2)*Utility::pow<2>(_polar_th_grad[_qp](1)) + 
                2.0*_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0)*_polar_th_grad[_qp](2) + 2.0*_grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1)*_polar_th_grad[_qp](2) + 3.0*_grad_phi[_j][_qp](2)*Utility::pow<2>(_polar_th_grad[_qp](2)) + 
                2.0*_alpha[_qp]*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](2) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1)*_polar_th_grad[_qp](2) + 2.0*_azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*
                 Utility::pow<3>((1.0/std::sin(_polar_th[_qp])))))*Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) + (std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*
           (-2.0*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*(_grad_test[_i][_qp](2) - _test[_i][_qp]) + 
             2.0*((_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2))*(_grad_test[_i][_qp](2) - _test[_i][_qp]) + 
                _alpha[_qp]*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2) - (_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_phi[_j][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*(1.0/std::sin(_polar_th[_qp])) + 
                _test[_i][_qp]*(-(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0)) - _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) - _grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2) + (Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_phi[_j][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*
                 Utility::pow<2>((1.0/std::sin(_polar_th[_qp]))))*Utility::pow<2>(std::sin(_azimuthal_ph[_qp]))) - _phi[_j][_qp]*Utility::pow<2>((1.0/std::sin(_polar_th[_qp])))*
           (-((Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*(_grad_test[_i][_qp](2) - _test[_i][_qp])) + 
             (_grad_test[_i][_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2))) - 
                (2.0*Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + 4.0*_azimuthal_ph_grad[_qp](1)*_azimuthal_ph_grad[_qp](2) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp] + 
                (1.0/std::sin(_polar_th[_qp]))*(2.0*_alpha[_qp]*_grad_test[_i][_qp](2)*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2)) - 
                   (Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2)))*_test[_i][_qp]*(1.0/std::sin(_polar_th[_qp]))))*Utility::pow<2>(std::sin(_azimuthal_ph[_qp])) + 
             (_azimuthal_ph_grad[_qp](0)*_grad_test[_i][_qp](0) + _azimuthal_ph_grad[_qp](2)*_grad_test[_i][_qp](1) + _azimuthal_ph_grad[_qp](1)*_grad_test[_i][_qp](2))*std::sin(2.0*_azimuthal_ph[_qp])) + 
          2.0*std::cos(_azimuthal_ph[_qp])*std::sin(_azimuthal_ph[_qp])*(2.0*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2))*_grad_test[_i][_qp](2) - 
             (3.0*_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*(2.0*_grad_phi[_j][_qp](1) + _grad_phi[_j][_qp](2)) + _azimuthal_ph_grad[_qp](2)*(_grad_phi[_j][_qp](1) + 2.0*_grad_phi[_j][_qp](2)))*_test[_i][_qp] + 
             2.0*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](2) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1)*_polar_th_grad[_qp](2) + 
                2.0*_azimuthal_ph_grad[_qp](2)*(_grad_phi[_j][_qp](0)*_polar_th_grad[_qp](0) + _grad_phi[_j][_qp](1)*_polar_th_grad[_qp](1) + 2.0*_grad_phi[_j][_qp](2)*_polar_th_grad[_qp](2)))*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])) - 
             2.0*(Utility::pow<3>(_azimuthal_ph_grad[_qp](2)) + (_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1))*_polar_th_grad[_qp](2) + 
                _azimuthal_ph_grad[_qp](2)*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + 2.0*Utility::pow<2>(_polar_th_grad[_qp](2))))*_phi[_j][_qp]*_test[_i][_qp]*Utility::pow<2>((1.0/std::sin(_polar_th[_qp]))) + 
             2.0*(_azimuthal_ph_grad[_qp](0)*_grad_phi[_j][_qp](0) + _azimuthal_ph_grad[_qp](1)*_grad_phi[_j][_qp](1) + _azimuthal_ph_grad[_qp](2)*_grad_phi[_j][_qp](2))*(-_grad_test[_i][_qp](2) + _test[_i][_qp] + _polar_th_grad[_qp](2)*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*Utility::pow<2>((1.0/std::sin(_polar_th[_qp]))) - 
             4.0*(_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_phi[_j][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp]))*(-_grad_test[_i][_qp](2) + _test[_i][_qp] + _polar_th_grad[_qp](2)*_test[_i][_qp]*(std::cos(_polar_th[_qp])/std::sin(_polar_th[_qp])))*Utility::pow<2>((1.0/std::sin(_polar_th[_qp]))) + 
             (_azimuthal_ph_grad[_qp](0)*_polar_th_grad[_qp](0) + _azimuthal_ph_grad[_qp](1)*_polar_th_grad[_qp](1) + _azimuthal_ph_grad[_qp](2)*_polar_th_grad[_qp](2))*_test[_i][_qp]*Utility::pow<4>((1.0/std::sin(_polar_th[_qp])))*(-2.0*_polar_th_grad[_qp](2)*_phi[_j][_qp] + _grad_phi[_j][_qp](2)*std::sin(2.0*_polar_th[_qp]))))))/((1.0 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
