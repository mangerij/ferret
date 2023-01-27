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

#include "AnisotropyCartLL.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", AnisotropyCartLL);

InputParameters AnisotropyCartLL::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  return params;
}

AnisotropyCartLL::AnisotropyCartLL(const InputParameters & parameters)
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
  _mu0(getMaterialProperty<Real>("mu0"))
{
}

Real
AnisotropyCartLL::computeQpResidual()
{
  if (_component == 0)
  {
    return -_mu0[_qp]*(-2*_g0[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp]))*(_mag_z[_qp]*_ny[_qp] - _mag_y[_qp]*_nz[_qp] + _alpha[_qp]*_Ms[_qp]*(-(Utility::pow<2>(_mag_y[_qp])*_nx[_qp]) + _mag_x[_qp]*_mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*(-(_mag_z[_qp]*_nx[_qp]) + _mag_x[_qp]*_nz[_qp])))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 1)
  {
    return -_mu0[_qp]*(2*_g0[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp]))*(_mag_z[_qp]*_nx[_qp] - _mag_x[_qp]*_nz[_qp] + _alpha[_qp]*_Ms[_qp]*(-(_mag_x[_qp]*_mag_y[_qp]*_nx[_qp]) + Utility::pow<2>(_mag_x[_qp])*_ny[_qp] + _mag_z[_qp]*(_mag_z[_qp]*_ny[_qp] - _mag_y[_qp]*_nz[_qp])))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 2)
  {
    return -_mu0[_qp]*(2*_g0[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(-(_mag_y[_qp]*_nx[_qp]) + _mag_x[_qp]*_ny[_qp] + _alpha[_qp]*_Ms[_qp]*(-(_mag_z[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp])) + (Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*_nz[_qp]))*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp]))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else
    return 0.0;
}

Real
AnisotropyCartLL::computeQpJacobian()
{
  if (_component == 0)
  {
    return -_mu0[_qp]*(2*_g0[_qp]*(-(_nx[_qp]*(_mag_z[_qp]*_ny[_qp] - _mag_y[_qp]*_nz[_qp])*(_K1[_qp] + 6*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp]))) + _alpha[_qp]*_Ms[_qp]*(-2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*
           (Utility::pow<2>(_mag_y[_qp])*(-3*Utility::pow<2>(_nx[_qp]) + Utility::pow<2>(_ny[_qp])) + 2*_mag_y[_qp]*_ny[_qp]*(2*_mag_x[_qp]*_nx[_qp] + _mag_z[_qp]*_nz[_qp]) + _mag_z[_qp]*(-3*_mag_z[_qp]*Utility::pow<2>(_nx[_qp]) + 4*_mag_x[_qp]*_nx[_qp]*_nz[_qp] + _mag_z[_qp]*Utility::pow<2>(_nz[_qp]))) + 
          _K1[_qp]*(Utility::pow<2>(_mag_y[_qp])*(_nx[_qp] - _ny[_qp])*(_nx[_qp] + _ny[_qp]) - 2*_mag_y[_qp]*_ny[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_z[_qp]*_nz[_qp]) + _mag_z[_qp]*(-2*_mag_x[_qp]*_nx[_qp]*_nz[_qp] + _mag_z[_qp]*(_nx[_qp] - _nz[_qp])*(_nx[_qp] + _nz[_qp])))))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 1)
  {
    return -_mu0[_qp]*(2*_g0[_qp]*(_ny[_qp]*(_mag_z[_qp]*_nx[_qp] - _mag_x[_qp]*_nz[_qp])*(_K1[_qp] + 6*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])) + _alpha[_qp]*_Ms[_qp]*(-(_K1[_qp]*(Utility::pow<2>(_mag_x[_qp])*(_nx[_qp] - _ny[_qp])*(_nx[_qp] + _ny[_qp]) + 2*_mag_x[_qp]*_nx[_qp]*(_mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp]) + _mag_z[_qp]*(-(_mag_z[_qp]*Utility::pow<2>(_ny[_qp])) + 2*_mag_y[_qp]*_ny[_qp]*_nz[_qp] + _mag_z[_qp]*Utility::pow<2>(_nz[_qp])))) - 
          2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(Utility::pow<2>(_mag_x[_qp])*(Utility::pow<2>(_nx[_qp]) - 3*Utility::pow<2>(_ny[_qp])) + 2*_mag_x[_qp]*_nx[_qp]*(2*_mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp]) + _mag_z[_qp]*(-3*_mag_z[_qp]*Utility::pow<2>(_ny[_qp]) + 4*_mag_y[_qp]*_ny[_qp]*_nz[_qp] + _mag_z[_qp]*Utility::pow<2>(_nz[_qp])))))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
  }
  else if (_component == 2)
  {
    return -_mu0[_qp]*(2*_g0[_qp]*(-(((_mag_y[_qp]*_nx[_qp] - _mag_x[_qp]*_ny[_qp])*_nz[_qp]*(_K1[_qp] + 6*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])))/_Ms[_qp]) + _alpha[_qp]*(-2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*
           (Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp]) + 4*_mag_z[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp])*_nz[_qp] - 3*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*Utility::pow<2>(_nz[_qp])) + _K1[_qp]*(-Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp]) - 2*_mag_z[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp])*_nz[_qp] + (Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*Utility::pow<2>(_nz[_qp]))))*_phi[_j][_qp]*_test[_i][_qp]
     )/(1 + Utility::pow<2>(_alpha[_qp]));
  }
  else
    return 0.0;
}

Real
AnisotropyCartLL::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return -_mu0[_qp]*(2*_g0[_qp]*((_alpha[_qp]*_Ms[_qp]*(2*_mag_y[_qp]*_nx[_qp] - _mag_x[_qp]*_ny[_qp]) + _nz[_qp])*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])) - 
       4*_K2[_qp]*_ny[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_mag_z[_qp]*_ny[_qp] - _mag_y[_qp]*_nz[_qp] + _alpha[_qp]*_Ms[_qp]*(-(Utility::pow<2>(_mag_y[_qp])*_nx[_qp]) + _mag_x[_qp]*_mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*(-(_mag_z[_qp]*_nx[_qp]) + _mag_x[_qp]*_nz[_qp]))) - 
       _ny[_qp]*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp]))*(_mag_z[_qp]*_ny[_qp] - _mag_y[_qp]*_nz[_qp] + _alpha[_qp]*_Ms[_qp]*(-(Utility::pow<2>(_mag_y[_qp])*_nx[_qp]) + _mag_x[_qp]*_mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*(-(_mag_z[_qp]*_nx[_qp]) + _mag_x[_qp]*_nz[_qp]))))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag_z_var)
    {
      return -_mu0[_qp]*(2*_g0[_qp]*(-((_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_ny[_qp] + _alpha[_qp]*_Ms[_qp]*(-2*_mag_z[_qp]*_nx[_qp] + _mag_x[_qp]*_nz[_qp]))*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp]))) - 
       4*_K2[_qp]*_nz[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_mag_z[_qp]*_ny[_qp] - _mag_y[_qp]*_nz[_qp] + _alpha[_qp]*_Ms[_qp]*(-(Utility::pow<2>(_mag_y[_qp])*_nx[_qp]) + _mag_x[_qp]*_mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*(-(_mag_z[_qp]*_nx[_qp]) + _mag_x[_qp]*_nz[_qp]))) - 
       _nz[_qp]*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp]))*(_mag_z[_qp]*_ny[_qp] - _mag_y[_qp]*_nz[_qp] + _alpha[_qp]*_Ms[_qp]*(-(Utility::pow<2>(_mag_y[_qp])*_nx[_qp]) + _mag_x[_qp]*_mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*(-(_mag_z[_qp]*_nx[_qp]) + _mag_x[_qp]*_nz[_qp]))))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
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
      return -_mu0[_qp]*(2*_g0[_qp]*(-((_alpha[_qp]*_Ms[_qp]*(_mag_y[_qp]*_nx[_qp] - 2*_mag_x[_qp]*_ny[_qp]) + _nz[_qp])*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp]))) + 
       4*_K2[_qp]*_nx[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_mag_z[_qp]*_nx[_qp] - _mag_x[_qp]*_nz[_qp] + _alpha[_qp]*_Ms[_qp]*(-(_mag_x[_qp]*_mag_y[_qp]*_nx[_qp]) + Utility::pow<2>(_mag_x[_qp])*_ny[_qp] + _mag_z[_qp]*(_mag_z[_qp]*_ny[_qp] - _mag_y[_qp]*_nz[_qp]))) + 
       _nx[_qp]*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp]))*(_mag_z[_qp]*_nx[_qp] - _mag_x[_qp]*_nz[_qp] + _alpha[_qp]*_Ms[_qp]*(-(_mag_x[_qp]*_mag_y[_qp]*_nx[_qp]) + Utility::pow<2>(_mag_x[_qp])*_ny[_qp] + _mag_z[_qp]*(_mag_z[_qp]*_ny[_qp] - _mag_y[_qp]*_nz[_qp]))))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag_z_var)
    {
      return -_mu0[_qp]*(2*_g0[_qp]*((_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_nx[_qp] + _alpha[_qp]*_Ms[_qp]*(2*_mag_z[_qp]*_ny[_qp] - _mag_y[_qp]*_nz[_qp]))*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])) + 
       4*_K2[_qp]*_nz[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_mag_z[_qp]*_nx[_qp] - _mag_x[_qp]*_nz[_qp] + _alpha[_qp]*_Ms[_qp]*(-(_mag_x[_qp]*_mag_y[_qp]*_nx[_qp]) + Utility::pow<2>(_mag_x[_qp])*_ny[_qp] + _mag_z[_qp]*(_mag_z[_qp]*_ny[_qp] - _mag_y[_qp]*_nz[_qp]))) + 
       _nz[_qp]*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp]))*(_mag_z[_qp]*_nx[_qp] - _mag_x[_qp]*_nz[_qp] + _alpha[_qp]*_Ms[_qp]*(-(_mag_x[_qp]*_mag_y[_qp]*_nx[_qp]) + Utility::pow<2>(_mag_x[_qp])*_ny[_qp] + _mag_z[_qp]*(_mag_z[_qp]*_ny[_qp] - _mag_y[_qp]*_nz[_qp]))))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
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
      return -_mu0[_qp]*(2*_g0[_qp]*(4*_K2[_qp]*_nx[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(-(_mag_y[_qp]*_nx[_qp]) + _mag_x[_qp]*_ny[_qp] + _alpha[_qp]*_Ms[_qp]*(-(_mag_z[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp])) + (Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*_nz[_qp])) + 
       (_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_ny[_qp] + _alpha[_qp]*_Ms[_qp]*(-(_mag_z[_qp]*_nx[_qp]) + 2*_mag_x[_qp]*_nz[_qp]))*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])) + 
       _nx[_qp]*(-(_mag_y[_qp]*_nx[_qp]) + _mag_x[_qp]*_ny[_qp] + _alpha[_qp]*_Ms[_qp]*(-(_mag_z[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp])) + (Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*_nz[_qp]))*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else if (jvar == _mag_y_var)
    {
      return -_mu0[_qp]*(2*_g0[_qp]*(4*_K2[_qp]*_ny[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(-(_mag_y[_qp]*_nx[_qp]) + _mag_x[_qp]*_ny[_qp] + _alpha[_qp]*_Ms[_qp]*(-(_mag_z[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp])) + (Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*_nz[_qp])) - 
       (_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])*(_nx[_qp] + _alpha[_qp]*_Ms[_qp]*(_mag_z[_qp]*_ny[_qp] - 2*_mag_y[_qp]*_nz[_qp]))*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])) + 
       _ny[_qp]*(-(_mag_y[_qp]*_nx[_qp]) + _mag_x[_qp]*_ny[_qp] + _alpha[_qp]*_Ms[_qp]*(-(_mag_z[_qp]*(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp])) + (Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*_nz[_qp]))*(_K1[_qp] + 2*_K2[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*_Ms[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
