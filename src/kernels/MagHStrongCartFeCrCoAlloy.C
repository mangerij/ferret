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

#include "MagHStrongCartFeCrCoAlloy.h"

class MagHStrongCartFeCrCoAlloy;

registerMooseObject("FerretApp", MagHStrongCartFeCrCoAlloy);

template<>
InputParameters validParams<MagHStrongCartFeCrCoAlloy>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for bound magnetic charge (div M)");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the magnetic vector");
  params.addRequiredCoupledVar("c1", "The concentration of Fe");
  params.addRequiredCoupledVar("c2", "The concentration of Cr");
  params.addRequiredCoupledVar("c3", "The concentration of Co");
  return params;
}

MagHStrongCartFeCrCoAlloy::MagHStrongCartFeCrCoAlloy(const InputParameters & parameters)
  :Kernel(parameters),
   _c2_var(coupled("c2")),
   _c3_var(coupled("c3")),
   _mag_x_var(coupled("mag_x")),
   _mag_y_var(coupled("mag_y")),
   _mag_z_var(coupled("mag_z")),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _c1(coupledValue("c1")),
   _c2(coupledValue("c2")),
   _c3(coupledValue("c3")),
   _bohrM(getMaterialProperty<Real>("bohrM")),
   _T(getMaterialProperty<Real>("T")),
   _mu0(getMaterialProperty<Real>("mu0")),
   _Ms(getMaterialProperty<Real>("Ms"))
{
}

Real
MagHStrongCartFeCrCoAlloy::computeQpResidual()
{
  return -_Ms[_qp]*_mu0[_qp]*(_grad_test[_i][_qp](0)*_mag_x[_qp]+_grad_test[_i][_qp](1)*_mag_y[_qp]+_grad_test[_i][_qp](2)*_mag_z[_qp]);
}
Real
MagHStrongCartFeCrCoAlloy::computeQpJacobian()
{
  return 0.0;
}

Real
MagHStrongCartFeCrCoAlloy::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _c2_var)
  {
    Real tau = _T[_qp]/(1043*_c1[_qp] - 311.5*_c2[_qp] + 1450*_c3[_qp] + (1650 + 550*(_c2[_qp] - _c1[_qp]))*_c1[_qp]*_c2[_qp] + 590*_c1[_qp]*_c3[_qp]);
    if (tau > 0.9)
    {
      return _mu0[_qp]*(_grad_test[_i][_qp](0)*_mag_x[_qp]+_grad_test[_i][_qp](1)*_mag_y[_qp]+_grad_test[_i][_qp](2)*_mag_z[_qp])*std::pow(2,8 + (0.01818181818181818*_T[_qp])/(0.5663636363636364*_c2[_qp] + _c1[_qp]*(-1.8963636363636365 + (-3. + _c1[_qp] - _c2[_qp])*_c2[_qp] - 1.0727272727272728*_c3[_qp]) - 2.6363636363636362*_c3[_qp])*_bohrM[_qp]*_phi[_j][_qp]*
   (-0.01 - 0.85*_c1[_qp] + (0.0030473270592617225*(0.5663636363636364 + _c1[_qp]*(-3. + _c1[_qp] - 2.*_c2[_qp]))*
        (0.0413564929693962*_c2[_qp] - 5.583126550868487*_c3[_qp] + _c1[_qp]*(-9.181141439205957 + 3.5153019023986767*_c2[_qp] + (-9.978081058726222 + _c1[_qp] - _c3[_qp])*_c3[_qp]))*_T[_qp])/
      Utility::pow<2>(0.5663636363636364*_c2[_qp] + _c1[_qp]*(-1.8963636363636365 + (-3. + _c1[_qp] - _c2[_qp])*_c2[_qp] - 1.0727272727272728*_c3[_qp]) - 2.6363636363636362*_c3[_qp])));
    }
    else if (tau <= 0.9)
    {
      return _mu0[_qp]*(_grad_test[_i][_qp](0)*_mag_x[_qp]+_grad_test[_i][_qp](1)*_mag_y[_qp]+_grad_test[_i][_qp](2)*_mag_z[_qp])*(20*_bohrM[_qp]*(-311.5 - 550.*Utility::pow<2>(_c1[_qp]) + _c1[_qp]*(1650. + 1100.*_c2[_qp]))*(2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*(2.4127 - 0.2418*_c1[_qp] + 0.2418*_c3[_qp])*_c3[_qp])*_phi[_j][_qp]*Utility::pow<4>(_T[_qp])*
      (Utility::pow<16>(1043.*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650. - 550.*_c1[_qp] + 550.*_c2[_qp]) + 1450.*_c3[_qp] + 590.*_c1[_qp]*_c3[_qp]) + 2*Utility::pow<16>(_T[_qp])))/(7.*Utility::pow<21>(1043.*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650. - 550.*_c1[_qp] + 550.*_c2[_qp]) + 1450.*_c3[_qp] + 590.*_c1[_qp]*_c3[_qp])) + 
   _bohrM[_qp]*(-0.01 - 0.85*_c1[_qp])*_phi[_j][_qp]*(1 - (5*Utility::pow<4>(_T[_qp]))/(7.*Utility::pow<4>(1043.*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650. - 550.*_c1[_qp] + 550.*_c2[_qp]) + 1450.*_c3[_qp] + 590.*_c1[_qp]*_c3[_qp])) - 
      (2*Utility::pow<20>(_T[_qp]))/(7.*Utility::pow<20>(1043.*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650. - 550.*_c1[_qp] + 550.*_c2[_qp]) + 1450.*_c3[_qp] + 590.*_c1[_qp]*_c3[_qp])));
    }
    else
      return 0.0;
  }
  else if (jvar == _c3_var)
  {
    Real tau = _T[_qp]/(1043*_c1[_qp] - 311.5*_c2[_qp] + 1450*_c3[_qp] + (1650 + 550*(_c2[_qp] - _c1[_qp]))*_c1[_qp]*_c2[_qp] + 590*_c1[_qp]*_c3[_qp]);
    if (tau > 0.9)
    {
      return _mu0[_qp]*(_grad_test[_i][_qp](0)*_mag_x[_qp]+_grad_test[_i][_qp](1)*_mag_y[_qp]+_grad_test[_i][_qp](2)*_mag_z[_qp])*std::pow(2,8 - (10*_T[_qp])/(1043.*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650. - 550.*_c1[_qp] + 550.*_c2[_qp]) + 1450.*_c3[_qp] + 590.*_c1[_qp]*_c3[_qp]))*_bohrM[_qp]*_phi[_j][_qp]*
   (1.35 + _c1[_qp]*(2.4127 - 0.2418*_c1[_qp] + 0.2418*_c3[_qp]) + 0.2418*_c1[_qp]*_c3[_qp] + (10*(1450. + 590.*_c1[_qp])*(2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*(2.4127 - 0.2418*_c1[_qp] + 0.2418*_c3[_qp])*_c3[_qp])*_T[_qp]*0.6931471805599453)/
      Utility::pow<2>(1043.*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650. - 550.*_c1[_qp] + 550.*_c2[_qp]) + 1450.*_c3[_qp] + 590.*_c1[_qp]*_c3[_qp]));
    }
    else if (tau <= 0.9)
    {
      return _mu0[_qp]*(_grad_test[_i][_qp](0)*_mag_x[_qp]+_grad_test[_i][_qp](1)*_mag_y[_qp]+_grad_test[_i][_qp](2)*_mag_z[_qp])*(20*_bohrM[_qp]*(1450. + 590.*_c1[_qp])*(2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*(2.4127 - 0.2418*_c1[_qp] + 0.2418*_c3[_qp])*_c3[_qp])*_phi[_j][_qp]*Utility::pow<4>(_T[_qp])*
      (Utility::pow<16>(1043.*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650. - 550.*_c1[_qp] + 550.*_c2[_qp]) + 1450.*_c3[_qp] + 590.*_c1[_qp]*_c3[_qp]) + 2*Utility::pow<16>(_T[_qp])))/(7.*Utility::pow<21>(1043.*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650. - 550.*_c1[_qp] + 550.*_c2[_qp]) + 1450.*_c3[_qp] + 590.*_c1[_qp]*_c3[_qp])) + 
   _bohrM[_qp]*(1.35 - 0.2418*Utility::pow<2>(_c1[_qp]) + _c1[_qp]*(2.4127 + 0.4836*_c3[_qp]))*_phi[_j][_qp]*(1 - (5*Utility::pow<4>(_T[_qp]))/(7.*Utility::pow<4>(1043.*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650. - 550.*_c1[_qp] + 550.*_c2[_qp]) + 1450.*_c3[_qp] + 590.*_c1[_qp]*_c3[_qp])) - 
      (2*Utility::pow<20>(_T[_qp]))/(7.*Utility::pow<20>(1043.*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650. - 550.*_c1[_qp] + 550.*_c2[_qp]) + 1450.*_c3[_qp] + 590.*_c1[_qp]*_c3[_qp])));
    }
    else
      return 0.0;
  } 
  else if (jvar == _mag_x_var)
  {
    return -_Ms[_qp]*_mu0[_qp]*(_grad_test[_i][_qp](0)*_phi[_j][_qp]);
  }
  else if (jvar == _mag_y_var)
  {
    return -_Ms[_qp]*_mu0[_qp]*(_grad_test[_i][_qp](1)*_phi[_j][_qp]);
  }
  else if (jvar == _mag_z_var)
  {
    return -_Ms[_qp]*_mu0[_qp]*(_grad_test[_i][_qp](2)*_phi[_j][_qp]);
  }
  else
    return 0.0;

}
