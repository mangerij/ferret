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

#include "LongitudinalLLB.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", LongitudinalLLB);


InputParameters LongitudinalLLB::validParams()

{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  params.addParam<Real>("g0", 1.0, "electron gyromagnetic factor");
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
  _alpha(getMaterialProperty<Real>("alpha")),
  _g0(getParam<Real>("g0")),
  _alpha_long(getMaterialProperty<Real>("alpha_long"))
{
}

Real
LongitudinalLLB::computeQpResidual()
{
   Real _temp1 = Utility::pow<2>(_mag_x[_qp])+Utility::pow<2>(_mag_y[_qp])+Utility::pow<2>(_mag_z[_qp]);
   Real _temp=_g0*_alpha_long[_qp]*_temp1*(_temp1-1.)/(1+Utility::pow<2>(_alpha[_qp]));
  if (_component == 0)
  {
   return 1.0*_temp*_mag_x[_qp]*_test[_i][_qp];
  }
  else if (_component == 1)
  {
   return 1.0*_temp*_mag_y[_qp]*_test[_i][_qp];
  }
  else if (_component == 2)
  {
   return 1.0*_temp*_mag_z[_qp]*_test[_i][_qp];
  }
  else
    return 0.0;
}

Real
LongitudinalLLB::computeQpJacobian()
{
  if (_component == 0)
  {
  return 1.0*_test[_i][_qp]*(_alpha_long[_qp]*_g0*(5*Utility::pow<4>(_mag_x[_qp]) + (-1 + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])) + Utility::pow<2>(_mag_x[_qp])*(-3 + 6*Utility::pow<2>(_mag_y[_qp]) + 6*Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha[_qp]));
  }
  else if (_component == 1)
  {
  return 1.0*_test[_i][_qp]*(_alpha_long[_qp]*_g0*(Utility::pow<4>(_mag_x[_qp]) + 5*Utility::pow<4>(_mag_y[_qp]) - Utility::pow<2>(_mag_z[_qp]) + Utility::pow<4>(_mag_z[_qp]) + Utility::pow<2>(_mag_x[_qp])*(-1 + 6*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp])) + Utility::pow<2>(_mag_y[_qp])*(-3 + 6*Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha[_qp]));
  }
  else if (_component == 2)
  {
  return 1.0*_test[_i][_qp]*(_alpha_long[_qp]*_g0*(Utility::pow<4>(_mag_x[_qp]) + Utility::pow<4>(_mag_y[_qp]) - 3*Utility::pow<2>(_mag_z[_qp]) + 5*Utility::pow<4>(_mag_z[_qp]) + Utility::pow<2>(_mag_y[_qp])*(-1 + 6*Utility::pow<2>(_mag_z[_qp])) + Utility::pow<2>(_mag_x[_qp])*(-1 + 2*Utility::pow<2>(_mag_y[_qp]) + 6*Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha[_qp]));
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
    return 1.0*_test[_i][_qp]*(2*_alpha_long[_qp]*_g0*_mag_x[_qp]*_mag_y[_qp]*(-1 + 2*Utility::pow<2>(_mag_x[_qp]) + 2*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
    return 1.0*_test[_i][_qp]*(2*_alpha_long[_qp]*_g0*_mag_x[_qp]*_mag_z[_qp]*(-1 + 2*Utility::pow<2>(_mag_x[_qp]) + 2*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha[_qp]));
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
    return 1.0*_test[_i][_qp]*(2*_alpha_long[_qp]*_g0*_mag_x[_qp]*_mag_y[_qp]*(-1 + 2*Utility::pow<2>(_mag_x[_qp]) + 2*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
    return 1.0*_test[_i][_qp]*(2*_alpha_long[_qp]*_g0*_mag_y[_qp]*_mag_z[_qp]*(-1 + 2*Utility::pow<2>(_mag_x[_qp]) + 2*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha[_qp]));
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
    return 1.0*_test[_i][_qp]*(2*_alpha_long[_qp]*_g0*_mag_x[_qp]*_mag_z[_qp]*(-1 + 2*Utility::pow<2>(_mag_x[_qp]) + 2*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (jvar == _mag_y_var)
    {
    return 1.0*_test[_i][_qp]*(2*_alpha_long[_qp]*_g0*_mag_y[_qp]*_mag_z[_qp]*(-1 + 2*Utility::pow<2>(_mag_x[_qp]) + 2*Utility::pow<2>(_mag_y[_qp]) + 2*Utility::pow<2>(_mag_z[_qp]))*_phi[_j][_qp])/(1 + Utility::pow<2>(_alpha[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
