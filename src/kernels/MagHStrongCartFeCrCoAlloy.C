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
  params.addRequiredCoupledVar("c1", "The concentration of Cr");
  params.addRequiredCoupledVar("c2", "The concentration of Cr");
  params.addRequiredCoupledVar("c3", "The concentration of Co");
  params.addRequiredParam<Real>("bohrM", "bohrM");
  params.addRequiredParam<Real>("T", "T");
  return params;
}

MagHStrongCartFeCrCoAlloy::MagHStrongCartFeCrCoAlloy(const InputParameters & parameters)
  :Kernel(parameters),
   _c2_var(coupled("c2")),
   _c3_var(coupled("c3")),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _c1(coupledValue("c1")),
   _c2(coupledValue("c2")),
   _c3(coupledValue("c3")),
   _bohrM(getParam<Real>("bohrM")),
   _T(getParam<Real>("T"))
{
}

Real
MagHStrongCartFeCrCoAlloy::computeQpResidual()
{
  return -(_grad_test[_i][_qp](0)*_mag_x[_qp]+_grad_test[_i][_qp](1)*_mag_y[_qp]+_grad_test[_i][_qp](2)*_mag_z[_qp]);
 // return -_Ms*(_grad_test[_i][_qp](0)*_mag_x[_qp]+_grad_test[_i][_qp](1)*_mag_y[_qp]+_grad_test[_i][_qp](2)*_mag_z[_qp]);
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
    Real tau = _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + 1450*_c3[_qp] + (1650 + 550*(_c2[_qp]-_c1[_qp]))*_c1[_qp]*_c2[_qp] + 590*_c1[_qp]*_c3[_qp]);
    if (tau > 0.9)
    {
      return -_grad_test[_i][_qp](0)*std::pow(2.0,-2 - 10*(-1 + _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))*_bohrM*(-0.01*_phi[_j][_qp] - 0.85*_c1[_qp]*_phi[_j][_qp]) + 
   (5*std::pow(2.0,-1.0 - 10*(-1 + _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))*(2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*_c3[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp])))*_bohrM*
      (-311.5*_phi[_j][_qp] + 550*_c1[_qp]*_c2[_qp]*_phi[_j][_qp] + _c1[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp]))*_phi[_j][_qp])*_T*std::log(2.0))/Utility::pow<2>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])-_grad_test[_i][_qp](1)*std::pow(2.0,-2 - 10*(-1 + _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))*_bohrM*(-0.01*_phi[_j][_qp] - 0.85*_c1[_qp]*_phi[_j][_qp]) + 
   (5*std::pow(2.0,-1.0 - 10*(-1 + _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))*(2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*_c3[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp])))*_bohrM*
      (-311.5*_phi[_j][_qp] + 550*_c1[_qp]*_c2[_qp]*_phi[_j][_qp] + _c1[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp]))*_phi[_j][_qp])*_T*std::log(2.0))/Utility::pow<2>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])-_grad_test[_i][_qp](2)*std::pow(2.0,-2 - 10*(-1 + _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))*_bohrM*(-0.01*_phi[_j][_qp] - 0.85*_c1[_qp]*_phi[_j][_qp]) + 
   (5*std::pow(2.0,-1.0 - 10*(-1 + _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))*(2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*_c3[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp])))*_bohrM*
      (-311.5*_phi[_j][_qp] + 550*_c1[_qp]*_c2[_qp]*_phi[_j][_qp] + _c1[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp]))*_phi[_j][_qp])*_T*std::log(2.0))/Utility::pow<2>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]);
    }
    else if (tau <= 0.9)
    {
      return -_grad_test[_i][_qp](0)*((2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*_c3[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp])))*_bohrM*((20*(-311.5*_phi[_j][_qp] + 550*_c1[_qp]*_c2[_qp]*_phi[_j][_qp] + _c1[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp]))*_phi[_j][_qp])*Utility::pow<4>(_T))/
         Utility::pow<5>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]) + 
        (40*(-311.5*_phi[_j][_qp] + 550*_c1[_qp]*_c2[_qp]*_phi[_j][_qp] + _c1[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp]))*_phi[_j][_qp])*Utility::pow<20>(_T))/Utility::pow<21>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))/7. + 
   _bohrM*(-0.01*_phi[_j][_qp] - 0.85*_c1[_qp]*_phi[_j][_qp])*(1 + ((-5*Utility::pow<4>(_T))/Utility::pow<4>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]) - 
         (2*Utility::pow<20>(_T))/Utility::pow<20>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]))/7.)-_grad_test[_i][_qp](1)*((2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*_c3[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp])))*_bohrM*((20*(-311.5*_phi[_j][_qp] + 550*_c1[_qp]*_c2[_qp]*_phi[_j][_qp] + _c1[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp]))*_phi[_j][_qp])*Utility::pow<4>(_T))/
         Utility::pow<5>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]) + 
        (40*(-311.5*_phi[_j][_qp] + 550*_c1[_qp]*_c2[_qp]*_phi[_j][_qp] + _c1[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp]))*_phi[_j][_qp])*Utility::pow<20>(_T))/Utility::pow<21>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))/7. + 
   _bohrM*(-0.01*_phi[_j][_qp] - 0.85*_c1[_qp]*_phi[_j][_qp])*(1 + ((-5*Utility::pow<4>(_T))/Utility::pow<4>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]) - 
         (2*Utility::pow<20>(_T))/Utility::pow<20>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]))/7.)-_grad_test[_i][_qp](2)*((2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*_c3[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp])))*_bohrM*((20*(-311.5*_phi[_j][_qp] + 550*_c1[_qp]*_c2[_qp]*_phi[_j][_qp] + _c1[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp]))*_phi[_j][_qp])*Utility::pow<4>(_T))/
         Utility::pow<5>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]) + 
        (40*(-311.5*_phi[_j][_qp] + 550*_c1[_qp]*_c2[_qp]*_phi[_j][_qp] + _c1[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp]))*_phi[_j][_qp])*Utility::pow<20>(_T))/Utility::pow<21>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))/7. + 
   _bohrM*(-0.01*_phi[_j][_qp] - 0.85*_c1[_qp]*_phi[_j][_qp])*(1 + ((-5*Utility::pow<4>(_T))/Utility::pow<4>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]) - 
         (2*Utility::pow<20>(_T))/Utility::pow<20>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]))/7.);
    }
    else
      return 0.0;
  }
  else if (jvar == _c3_var)
  {
    Real tau = _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + 1450*_c3[_qp] + (1650 + 550*(_c2[_qp]-_c1[_qp]))*_c1[_qp]*_c2[_qp] + 590*_c1[_qp]*_c3[_qp]);
    if (tau > 0.9)
    {
      return -_grad_test[_i][_qp](0)*std::pow(2.0,-2 - 10*(-1 + _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))*_bohrM*(1.35*_phi[_j][_qp] + 0.2418*_c1[_qp]*_c3[_qp]*_phi[_j][_qp] + _c1[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp]))*_phi[_j][_qp]) + 
   (5*std::pow(2.0,-1.0 - 10*(-1 + _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))*(2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*_c3[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp])))*_bohrM*(1450*_phi[_j][_qp] + 590*_c1[_qp]*_phi[_j][_qp])*
      _T*std::log(2.0))/Utility::pow<2>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])-_grad_test[_i][_qp](1)*std::pow(2.0,-2 - 10*(-1 + _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))*_bohrM*(1.35*_phi[_j][_qp] + 0.2418*_c1[_qp]*_c3[_qp]*_phi[_j][_qp] + _c1[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp]))*_phi[_j][_qp]) + 
   (5*std::pow(2.0,-1.0 - 10*(-1 + _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))*(2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*_c3[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp])))*_bohrM*(1450*_phi[_j][_qp] + 590*_c1[_qp]*_phi[_j][_qp])*
      _T*std::log(2.0))/Utility::pow<2>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])-_grad_test[_i][_qp](2)*std::pow(2.0,-2 - 10*(-1 + _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))*_bohrM*(1.35*_phi[_j][_qp] + 0.2418*_c1[_qp]*_c3[_qp]*_phi[_j][_qp] + _c1[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp]))*_phi[_j][_qp]) + 
   (5*std::pow(2.0,-1.0 - 10*(-1 + _T/(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))*(2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*_c3[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp])))*_bohrM*(1450*_phi[_j][_qp] + 590*_c1[_qp]*_phi[_j][_qp])*
      _T*std::log(2.0))/Utility::pow<2>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]);
    }
    else if (tau <= 0.9)
    {
      return -_grad_test[_i][_qp](0)*((2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*_c3[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp])))*_bohrM*((20*(1450*_phi[_j][_qp] + 590*_c1[_qp]*_phi[_j][_qp])*Utility::pow<4>(_T))/Utility::pow<5>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]) + 
        (40*(1450*_phi[_j][_qp] + 590*_c1[_qp]*_phi[_j][_qp])*Utility::pow<20>(_T))/Utility::pow<21>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))/7. + 
   _bohrM*(1.35*_phi[_j][_qp] + 0.2418*_c1[_qp]*_c3[_qp]*_phi[_j][_qp] + _c1[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp]))*_phi[_j][_qp])*(1 + ((-5*Utility::pow<4>(_T))/Utility::pow<4>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]) - 
         (2*Utility::pow<20>(_T))/Utility::pow<20>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]))/7.)-_grad_test[_i][_qp](1)*((2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*_c3[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp])))*_bohrM*((20*(1450*_phi[_j][_qp] + 590*_c1[_qp]*_phi[_j][_qp])*Utility::pow<4>(_T))/Utility::pow<5>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]) + 
        (40*(1450*_phi[_j][_qp] + 590*_c1[_qp]*_phi[_j][_qp])*Utility::pow<20>(_T))/Utility::pow<21>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))/7. + 
   _bohrM*(1.35*_phi[_j][_qp] + 0.2418*_c1[_qp]*_c3[_qp]*_phi[_j][_qp] + _c1[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp]))*_phi[_j][_qp])*(1 + ((-5*Utility::pow<4>(_T))/Utility::pow<4>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]) - 
         (2*Utility::pow<20>(_T))/Utility::pow<20>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]))/7.)-_grad_test[_i][_qp](2)*((2.22*_c1[_qp] - 0.01*_c2[_qp] - 0.85*_c1[_qp]*_c2[_qp] + 1.35*_c3[_qp] + _c1[_qp]*_c3[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp])))*_bohrM*((20*(1450*_phi[_j][_qp] + 590*_c1[_qp]*_phi[_j][_qp])*Utility::pow<4>(_T))/Utility::pow<5>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]) + 
        (40*(1450*_phi[_j][_qp] + 590*_c1[_qp]*_phi[_j][_qp])*Utility::pow<20>(_T))/Utility::pow<21>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp])))/7. + 
   _bohrM*(1.35*_phi[_j][_qp] + 0.2418*_c1[_qp]*_c3[_qp]*_phi[_j][_qp] + _c1[_qp]*(2.4127 + 0.2418*(-_c1[_qp] + _c3[_qp]))*_phi[_j][_qp])*(1 + ((-5*Utility::pow<4>(_T))/Utility::pow<4>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]) - 
         (2*Utility::pow<20>(_T))/Utility::pow<20>(1043*_c1[_qp] - 311.5*_c2[_qp] + _c1[_qp]*_c2[_qp]*(1650 + 550*(-_c1[_qp] + _c2[_qp])) + 1450*_c3[_qp] + 590*_c1[_qp]*_c3[_qp]))/7.);
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}
