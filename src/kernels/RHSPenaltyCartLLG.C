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

#include "RHSPenaltyCartLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", RHSPenaltyCartLLG);

template<>
InputParameters validParams<RHSPenaltyCartLLG>()
{
  InputParameters params = validParams<TimeKernel>();
  params.addClassDescription("Calculates a residual contribution - M$*$H in the total energy, assuming H = - div * potential.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  params.addRequiredParam<Real>("g0", "g0");
  params.addRequiredParam<Real>("mu0", "mu0");
  params.addRequiredParam<Real>("Ap", "Ap");
  params.addRequiredParam<Real>("eps", "eps");
  return params;
}

RHSPenaltyCartLLG::RHSPenaltyCartLLG(const InputParameters & parameters)
  :TimeKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _mag_x_var(coupled("mag_x")),
   _mag_y_var(coupled("mag_y")),
   _mag_z_var(coupled("mag_z")),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _mag_x_dot(coupledDot("mag_x")),
   _mag_y_dot(coupledDot("mag_y")),
   _mag_z_dot(coupledDot("mag_z")),
   _mag_x_d_dot(coupledDotDu("mag_x")),
   _mag_y_d_dot(coupledDotDu("mag_y")),
   _mag_z_d_dot(coupledDotDu("mag_z")),
   _g0(getParam<Real>("g0")),
   _mu0(getParam<Real>("mu0")),
   _Ap(getParam<Real>("Ap")),
   _eps(getParam<Real>("eps"))
{
}


Real
RHSPenaltyCartLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp]*_mu0*_g0*_Ap*(_mag_y[_qp]*_mag_z_dot[_qp]-_mag_z[_qp]*_mag_y_dot[_qp]+_eps*_mag_x[_qp]);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp]*_mu0*_g0*_Ap*(_mag_z[_qp]*_mag_x_dot[_qp]-_mag_x[_qp]*_mag_z_dot[_qp]+_eps*_mag_y[_qp]);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp]*_mu0*_g0*_Ap*(_mag_x[_qp]*_mag_y_dot[_qp]-_mag_y[_qp]*_mag_x_dot[_qp]+_eps*_mag_z[_qp]);
  }
  else
    return 0.0;
}

Real
RHSPenaltyCartLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp]*_phi[_j][_qp]*_eps*_mu0*_g0*_Ap;
  }
  else if (_component == 1)
  {
    return _test[_i][_qp]*_phi[_j][_qp]*_eps*_mu0*_g0*_Ap;
  }
  else if (_component == 2)
  {
    return _test[_i][_qp]*_phi[_j][_qp]*_eps*_mu0*_g0*_Ap;
  }
  else
    return 0.0;
}

Real
RHSPenaltyCartLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return _test[_i][_qp]*_phi[_j][_qp]*(_mag_z_dot[_qp]-_mag_z[_qp]*_mag_y_d_dot[_qp]);
    }
    else if (jvar == _mag_z_var)
    {
      return _test[_i][_qp]*_phi[_j][_qp]*(-_mag_y_dot[_qp]+_mag_y[_qp]*_mag_z_d_dot[_qp]);
    }
    else
      return 0.0;
  }
  else if (_component == 1)
  {
    if (jvar == _mag_x_var)
    {
      return _test[_i][_qp]*_phi[_j][_qp]*(-_mag_z_dot[_qp]+_mag_z[_qp]*_mag_x_d_dot[_qp]);
    }
    else if (jvar == _mag_z_var)
    {
      return _test[_i][_qp]*_phi[_j][_qp]*(_mag_x_dot[_qp]-_mag_x[_qp]*_mag_z_d_dot[_qp]);
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
      return _test[_i][_qp]*_phi[_j][_qp]*(_mag_y_dot[_qp]-_mag_y[_qp]*_mag_x_d_dot[_qp]);
    }
    else if (jvar == _mag_y_var)
    {
      return _test[_i][_qp]*_phi[_j][_qp]*(-_mag_x_dot[_qp]+_mag_x[_qp]*_mag_y_d_dot[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
