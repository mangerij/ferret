/*
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the ter___Ms of the GNU General Public License as published by
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

#include "LongitudinalLLB.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", LongitudinalLLB);

template<>
InputParameters validParams<LongitudinalLLB>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  params.addRequiredParam<Real>("alpha", "the damping coefficient in the LLG equation");
  params.addRequiredParam<Real>("g0", "g0");
  params.addRequiredParam<Real>("Ms", "Ms");
  params.addRequiredParam<Real>("Ae", "Ae");
  params.addRequiredParam<Real>("mu0", "mu0");
  params.addRequiredParam<Real>("alpha_long", "alpha_long");
  return params;
}

LongitudinalLLB::LongitudinalLLB(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _alpha(getParam<Real>("alpha")),
  _g0(getParam<Real>("g0")),
  _Ae(getParam<Real>("Ae")),
  _Ms(getParam<Real>("Ms")),
  _mu0(getParam<Real>("mu0")),
  _alpha_long(getParam<Real>("alpha_long"))
{
}

Real
LongitudinalLLB::computeQpResidual()
{
   Real _temp1 = Utility::pow<2>(_mag_x[_qp])+Utility::pow<2>(_mag_y[_qp])+Utility::pow<2>(_mag_z[_qp]);
   Real _temp=_g0*_alpha_long*_temp1*(_temp1-1.)/(1+Utility::pow<2>(_alpha));
  if (_component == 0)
  {
   return _temp*_mag_x[_qp]*_test[_i][_qp];
  }
  else if (_component == 1)
  {
   return _temp*_mag_y[_qp]*_test[_i][_qp];
  }
  else if (_component == 2)
  {
   return _temp*_mag_z[_qp]*_test[_i][_qp];
  }
  else
    return 0.0;
}

Real
LongitudinalLLB::computeQpJacobian()
{
  if (_component == 0)
  {
  return  _test[_i][_qp]*(_alpha_long*_g0*(5*Utility::pow<4>(_mag_x[_qp]) + (-1 + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])) + Utility::pow<2>(_mag_x[_qp])*(-3 + 6*Utility::pow<2>(_mag_y[_qp]) + 6*Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha));
  }
  else if (_component == 1)
  {
  return _test[_i][_qp]*(_alpha_long*_g0*(Utility::pow<4>(_mag_x[_qp]) + 5*Utility::pow<4>(_mag_y[_qp]) - Utility::pow<2>(_mag_z[_qp]) + Utility::pow<4>(_mag_z[_qp]) + Utility::pow<2>(_mag_x[_qp])*(-1 + 6*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp])) + Utility::pow<2>(_mag_y[_qp])*(-3 + 6*Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha));
  }
  else if (_component == 2)
  {
  return _test[_i][_qp]*(_alpha_long*_g0*(Utility::pow<4>(_mag_x[_qp]) + Utility::pow<4>(_mag_y[_qp]) - 3*Utility::pow<2>(_mag_z[_qp]) + 5*Utility::pow<4>(_mag_z[_qp]) + Utility::pow<2>(_mag_y[_qp])*(-1 + 6*Utility::pow<2>(_mag_z[_qp])) + Utility::pow<2>(_mag_x[_qp])*(-1 + 2*Utility::pow<2>(_mag_y[_qp]) + 6*Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha));
  }
  else
    return 0.0;
}

Real
LongitudinalLLB::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
    return _test[_i][_qp]*(2*_alpha_long*_g0*_mag_x[_qp]*_mag_y[_qp]*(-1 + 2*Utility::pow<2>(_mag_x[_qp]) + 2*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha));
    }
    else if (jvar == _mag_z_var)
    {
    return _test[_i][_qp]*(2*_alpha_long*_g0*_mag_x[_qp]*_mag_z[_qp]*(-1 + 2*Utility::pow<2>(_mag_x[_qp]) + 2*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _mag_x_var)
    {
    return _test[_i][_qp]*(2*_alpha_long*_g0*_mag_x[_qp]*_mag_y[_qp]*(-1 + 2*Utility::pow<2>(_mag_x[_qp]) + 2*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha)); 
    }
    else if (jvar == _mag_z_var)
    {
    return _test[_i][_qp]*(2*_alpha_long*_g0*_mag_y[_qp]*_mag_z[_qp]*(-1 + 2*Utility::pow<2>(_mag_x[_qp]) + 2*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _mag_x_var)
    {
    return _test[_i][_qp]*(2*_alpha_long*_g0*_mag_x[_qp]*_mag_z[_qp]*(-1 + 2*Utility::pow<2>(_mag_x[_qp]) + 2*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha));
    }
    else if (jvar == _mag_y_var)
    {
    return _test[_i][_qp]*(2*_alpha_long*_g0*_mag_y[_qp]*_mag_z[_qp]*(-1 + 2*Utility::pow<2>(_mag_x[_qp]) + 2*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
