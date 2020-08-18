/*
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) a_ny[_qp] later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT A_ny[_qp] WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "LongitudinalAnisotropyCartLLB.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", LongitudinalAnisotropyCartLLB);

template<>
InputParameters validParams<LongitudinalAnisotropyCartLLB>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  return params;
}

LongitudinalAnisotropyCartLLB::LongitudinalAnisotropyCartLLB(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _mag_x_var(coupled("mag_x")),
  _mag_y_var(coupled("mag_y")),
  _mag_z_var(coupled("mag_z")),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _mag_x_grad(coupledGradient("mag_x")),
  _mag_y_grad(coupledGradient("mag_y")),
  _mag_z_grad(coupledGradient("mag_z")),
  _alpha(getMaterialProperty<Real>("alpha")),
  _K1(getMaterialProperty<Real>("K1")),
  _K2(getMaterialProperty<Real>("K2")),
  _nx(getMaterialProperty<Real>("nx")),
  _ny(getMaterialProperty<Real>("ny")),
  _nz(getMaterialProperty<Real>("nz")),
  _g0(getMaterialProperty<Real>("g0")),
  _Ms(getMaterialProperty<Real>("Ms")),
  _alpha_long(getMaterialProperty<Real>("alpha_long"))
{
}

Real
LongitudinalAnisotropyCartLLB::computeQpResidual()
{
  if (_component == 0)
  {
    return (-2*_alpha_long[_qp]*_g0[_qp]*_K1[_qp]*_mag_x[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<4>(_Ms[_qp]));
  }
  else if (_component == 1)
  {
    return (-2*_alpha_long[_qp]*_g0[_qp]*_K1[_qp]*_mag_y[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<4>(_Ms[_qp]));
  }
  else if (_component == 2)
  {
    return (-2*_alpha_long[_qp]*_g0[_qp]*_K1[_qp]*_mag_z[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<4>(_Ms[_qp]));
  }
  else
    return 0.0;
}

Real
LongitudinalAnisotropyCartLLB::computeQpJacobian()
{
  if (_component == 0)
  {
    return (-2*_alpha_long[_qp]*_g0[_qp]*_K1[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(3*_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<4>(_Ms[_qp]));
  }
  else if (_component == 1)
  {
    return (-2*_alpha_long[_qp]*_g0[_qp]*_K1[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_mag_x[_qp]*_nx[_qp] + 3*_mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<4>(_Ms[_qp]));
  }
  else if (_component == 2)
  {
    return (-2*_alpha_long[_qp]*_g0[_qp]*_K1[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + 3*_mag_z[_qp]*_nz[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<4>(_Ms[_qp]));
  }
  else
    return 0.0;
}

Real
LongitudinalAnisotropyCartLLB::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return (-4*_alpha_long[_qp]*_g0[_qp]*_K1[_qp]*_mag_x[_qp]*_ny[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<4>(_Ms[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
      return (-4*_alpha_long[_qp]*_g0[_qp]*_K1[_qp]*_mag_x[_qp]*_nz[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<4>(_Ms[_qp]));
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
      return (-4*_alpha_long[_qp]*_g0[_qp]*_K1[_qp]*_mag_y[_qp]*_nx[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<4>(_Ms[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
      return (-4*_alpha_long[_qp]*_g0[_qp]*_K1[_qp]*_mag_y[_qp]*_nz[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<4>(_Ms[_qp]));
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
      return (-4*_alpha_long[_qp]*_g0[_qp]*_K1[_qp]*_mag_z[_qp]*_nx[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<4>(_Ms[_qp]));
    }
    else if (jvar == _mag_y_var)
    {
      return (-4*_alpha_long[_qp]*_g0[_qp]*_K1[_qp]*_mag_z[_qp]*_ny[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<4>(_Ms[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
