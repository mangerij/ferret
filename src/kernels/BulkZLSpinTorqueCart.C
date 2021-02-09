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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "BulkZLSpinTorqueCart.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", BulkZLSpinTorqueCart);

template<>
InputParameters validParams<BulkZLSpinTorqueCart>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  return params;
}

BulkZLSpinTorqueCart::BulkZLSpinTorqueCart(const InputParameters & parameters)
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
  _j_x(coupledValue("j_x")),
  _j_y(coupledValue("j_y")),
  _j_z(coupledValue("j_z")),
  _muB(getMaterialProperty<Real>("muB")),
  _g0(getMaterialProperty<Real>("g0")),
  _e(getMaterialProperty<Real>("e")),
  _alpha(getMaterialProperty<Real>("alpha")),
  _Ms(getMaterialProperty<Real>("Ms")),
  _xi(getMaterialProperty<Real>("xi"))
{
}

Real
BulkZLSpinTorqueCart::computeQpResidual()
{
  if (_component == 0)
  {
    return (_muB[_qp]*_test[_i][_qp]*(_mag_y_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp]*_mag_y[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp]*_mag_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp]*_mag_y[_qp] - _mag_x_grad[_qp](0)*_j_x[_qp]*Utility::pow<2>(_mag_y[_qp]) - _mag_x_grad[_qp](1)*_j_y[_qp]*Utility::pow<2>(_mag_y[_qp]) - _mag_x_grad[_qp](2)*_j_z[_qp]*Utility::pow<2>(_mag_y[_qp]) + _mag_z_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp]*_mag_z[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp]*_mag_z[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp]*_mag_z[_qp] - _mag_x_grad[_qp](0)*_j_x[_qp]*Utility::pow<2>(_mag_z[_qp]) - 
       _mag_x_grad[_qp](1)*_j_y[_qp]*Utility::pow<2>(_mag_z[_qp]) - _mag_x_grad[_qp](2)*_j_z[_qp]*Utility::pow<2>(_mag_z[_qp]) + ((_mag_z_grad[_qp](0)*_j_x[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp])*_mag_y[_qp] - (_mag_y_grad[_qp](0)*_j_x[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp])*_mag_z[_qp])*_xi[_qp] + 
       _alpha[_qp]*(-(_mag_z_grad[_qp](0)*_j_x[_qp]*_mag_y[_qp]) - _mag_z_grad[_qp](1)*_j_y[_qp]*_mag_y[_qp] - _mag_z_grad[_qp](2)*_j_z[_qp]*_mag_y[_qp] + _mag_y_grad[_qp](0)*_j_x[_qp]*_mag_z[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp]*_mag_z[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp]*_mag_z[_qp] + 
          (_mag_y[_qp]*((_mag_y_grad[_qp](0)*_j_x[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp])*_mag_x[_qp] - (_mag_x_grad[_qp](0)*_j_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp])*_mag_y[_qp]) + (_mag_z_grad[_qp](0)*_j_x[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp])*_mag_x[_qp]*_mag_z[_qp] - (_mag_x_grad[_qp](0)*_j_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp])*Utility::pow<2>(_mag_z[_qp]))*_xi[_qp])))/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*_e[_qp]*_g0[_qp]*_Ms[_qp]*(1 + Utility::pow<2>(_xi[_qp])));
  }
  else if (_component == 1)
  {
    return (_muB[_qp]*_test[_i][_qp]*(-(_mag_y_grad[_qp](1)*_j_y[_qp]*Utility::pow<2>(_mag_x[_qp])) - _mag_y_grad[_qp](2)*_j_z[_qp]*Utility::pow<2>(_mag_x[_qp]) + _mag_x_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp]*_mag_y[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp]*_mag_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp]*_mag_y[_qp] + _mag_z_grad[_qp](0)*_j_x[_qp]*_mag_y[_qp]*_mag_z[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp]*_mag_y[_qp]*_mag_z[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp]*_mag_y[_qp]*_mag_z[_qp] - _mag_y_grad[_qp](1)*_j_y[_qp]*Utility::pow<2>(_mag_z[_qp]) - _mag_y_grad[_qp](2)*_j_z[_qp]*Utility::pow<2>(_mag_z[_qp]) - 
       _mag_y_grad[_qp](0)*_j_x[_qp]*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_z[_qp])) + (-((_mag_z_grad[_qp](0)*_j_x[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp])*_mag_x[_qp]) + (_mag_x_grad[_qp](0)*_j_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp])*_mag_z[_qp])*_xi[_qp] + 
       _alpha[_qp]*(_mag_z_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp] - _mag_x_grad[_qp](0)*_j_x[_qp]*_mag_z[_qp] - _mag_x_grad[_qp](1)*_j_y[_qp]*_mag_z[_qp] - _mag_x_grad[_qp](2)*_j_z[_qp]*_mag_z[_qp] + 
          (_mag_x[_qp]*(-((_mag_y_grad[_qp](0)*_j_x[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp])*_mag_x[_qp]) + (_mag_x_grad[_qp](0)*_j_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp])*_mag_y[_qp]) + (_mag_z_grad[_qp](0)*_j_x[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp])*_mag_y[_qp]*_mag_z[_qp] - (_mag_y_grad[_qp](0)*_j_x[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp])*Utility::pow<2>(_mag_z[_qp]))*_xi[_qp])))/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*_e[_qp]*_g0[_qp]*_Ms[_qp]*(1 + Utility::pow<2>(_xi[_qp])));
  }
  else if (_component == 2)
  {
    return (_muB[_qp]*_test[_i][_qp]*(-(_mag_z_grad[_qp](1)*_j_y[_qp]*Utility::pow<2>(_mag_x[_qp])) - _mag_z_grad[_qp](2)*_j_z[_qp]*Utility::pow<2>(_mag_x[_qp]) - _mag_z_grad[_qp](1)*_j_y[_qp]*Utility::pow<2>(_mag_y[_qp]) - _mag_z_grad[_qp](2)*_j_z[_qp]*Utility::pow<2>(_mag_y[_qp]) - _mag_z_grad[_qp](0)*_j_x[_qp]*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp])) + _mag_x_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp]*_mag_z[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp]*_mag_z[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp]*_mag_z[_qp] + _mag_y_grad[_qp](0)*_j_x[_qp]*_mag_y[_qp]*_mag_z[_qp] + 
       _mag_y_grad[_qp](1)*_j_y[_qp]*_mag_y[_qp]*_mag_z[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp]*_mag_y[_qp]*_mag_z[_qp] + ((_mag_y_grad[_qp](0)*_j_x[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp])*_mag_x[_qp] - (_mag_x_grad[_qp](0)*_j_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp])*_mag_y[_qp])*_xi[_qp] + 
       _alpha[_qp]*(-(_mag_y_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp]) - _mag_y_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp] - _mag_y_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp] + _mag_x_grad[_qp](0)*_j_x[_qp]*_mag_y[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp]*_mag_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp]*_mag_y[_qp] + 
          (-((_mag_z_grad[_qp](0)*_j_x[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp])*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))) + (_mag_x_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp] + _mag_y_grad[_qp](0)*_j_x[_qp]*_mag_y[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp]*_mag_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp]*_mag_y[_qp])*_mag_z[_qp])*_xi[_qp])))/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*_e[_qp]*_g0[_qp]*_Ms[_qp]*(1 + Utility::pow<2>(_xi[_qp])));
  }
  else
    return 0.0;
}

Real
BulkZLSpinTorqueCart::computeQpJacobian()
{
  if (_component == 0)
  {
    return -(_muB[_qp]*((_grad_phi[_j][_qp](0)*_j_x[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp])*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])) - (_mag_y_grad[_qp](0)*_j_x[_qp]*_mag_y[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp]*_mag_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp]*_mag_y[_qp] + _mag_z_grad[_qp](0)*_j_x[_qp]*_mag_z[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp]*_mag_z[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp]*_mag_z[_qp])*_phi[_j][_qp])*_test[_i][_qp]*(1 + _alpha[_qp]*_xi[_qp]))/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*_e[_qp]*_g0[_qp]*_Ms[_qp]*(1 + Utility::pow<2>(_xi[_qp])));
  }
  else if (_component == 1)
  {
    return -(_muB[_qp]*((_grad_phi[_j][_qp](0)*_j_x[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp])*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_z[_qp])) - (_mag_x_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp] + _mag_z_grad[_qp](0)*_j_x[_qp]*_mag_z[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp]*_mag_z[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp]*_mag_z[_qp])*_phi[_j][_qp])*_test[_i][_qp]*(1 + _alpha[_qp]*_xi[_qp]))/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*_e[_qp]*_g0[_qp]*_Ms[_qp]*(1 + Utility::pow<2>(_xi[_qp])));
  }
  else if (_component == 2)
  {
    return -(_muB[_qp]*((_grad_phi[_j][_qp](0)*_j_x[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp])*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp])) - (_mag_x_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp] + _mag_y_grad[_qp](0)*_j_x[_qp]*_mag_y[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp]*_mag_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp]*_mag_y[_qp])*_phi[_j][_qp])*_test[_i][_qp]*(1 + _alpha[_qp]*_xi[_qp]))/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*_e[_qp]*_g0[_qp]*_Ms[_qp]*(1 + Utility::pow<2>(_xi[_qp])));
  }
  else
    return 0.0;
}

Real
BulkZLSpinTorqueCart::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return (_muB[_qp]*_test[_i][_qp]*(_grad_phi[_j][_qp](0)*_j_x[_qp]*_mag_x[_qp]*_mag_y[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp]*_mag_x[_qp]*_mag_y[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp]*_mag_x[_qp]*_mag_y[_qp] + _alpha[_qp]*_grad_phi[_j][_qp](0)*_j_x[_qp]*_mag_z[_qp] + _alpha[_qp]*_grad_phi[_j][_qp](1)*_j_y[_qp]*_mag_z[_qp] + _alpha[_qp]*_grad_phi[_j][_qp](2)*_j_z[_qp]*_mag_z[_qp] - _alpha[_qp]*_mag_z_grad[_qp](0)*_j_x[_qp]*_phi[_j][_qp] - _alpha[_qp]*_mag_z_grad[_qp](1)*_j_y[_qp]*_phi[_j][_qp] - _alpha[_qp]*_mag_z_grad[_qp](2)*_j_z[_qp]*_phi[_j][_qp] + _mag_y_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp]*_phi[_j][_qp] + 
       _mag_y_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp]*_phi[_j][_qp] + _mag_y_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp]*_phi[_j][_qp] - 2*_mag_x_grad[_qp](0)*_j_x[_qp]*_mag_y[_qp]*_phi[_j][_qp] - 2*_mag_x_grad[_qp](1)*_j_y[_qp]*_mag_y[_qp]*_phi[_j][_qp] - 2*_mag_x_grad[_qp](2)*_j_z[_qp]*_mag_y[_qp]*_phi[_j][_qp] + 
       ((_grad_phi[_j][_qp](0)*_j_x[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp])*(_alpha[_qp]*_mag_x[_qp]*_mag_y[_qp] - _mag_z[_qp]) + (_mag_z_grad[_qp](0)*_j_x[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp] + _alpha[_qp]*(_mag_y_grad[_qp](0)*_j_x[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp])*_mag_x[_qp] - 2*_alpha[_qp]*(_mag_x_grad[_qp](0)*_j_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp])*_mag_y[_qp])*_phi[_j][_qp])*_xi[_qp]))/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*_e[_qp]*_g0[_qp]*_Ms[_qp]*(1 + Utility::pow<2>(_xi[_qp])));
    }
    else if (jvar == _mag_z_var)
    {
      return (_muB[_qp]*_test[_i][_qp]*(-(_alpha[_qp]*_grad_phi[_j][_qp](0)*_j_x[_qp]*_mag_y[_qp]) - _alpha[_qp]*_grad_phi[_j][_qp](1)*_j_y[_qp]*_mag_y[_qp] - _alpha[_qp]*_grad_phi[_j][_qp](2)*_j_z[_qp]*_mag_y[_qp] + _grad_phi[_j][_qp](0)*_j_x[_qp]*_mag_x[_qp]*_mag_z[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp]*_mag_x[_qp]*_mag_z[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp]*_mag_x[_qp]*_mag_z[_qp] + _alpha[_qp]*_mag_y_grad[_qp](0)*_j_x[_qp]*_phi[_j][_qp] + _alpha[_qp]*_mag_y_grad[_qp](1)*_j_y[_qp]*_phi[_j][_qp] + _alpha[_qp]*_mag_y_grad[_qp](2)*_j_z[_qp]*_phi[_j][_qp] + _mag_z_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp]*_phi[_j][_qp] + 
       _mag_z_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp]*_phi[_j][_qp] + _mag_z_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp]*_phi[_j][_qp] - 2*_mag_x_grad[_qp](0)*_j_x[_qp]*_mag_z[_qp]*_phi[_j][_qp] - 2*_mag_x_grad[_qp](1)*_j_y[_qp]*_mag_z[_qp]*_phi[_j][_qp] - 2*_mag_x_grad[_qp](2)*_j_z[_qp]*_mag_z[_qp]*_phi[_j][_qp] + 
       ((_grad_phi[_j][_qp](0)*_j_x[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp])*(_mag_y[_qp] + _alpha[_qp]*_mag_x[_qp]*_mag_z[_qp]) - (_mag_y_grad[_qp](0)*_j_x[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp] - _alpha[_qp]*(_mag_z_grad[_qp](0)*_j_x[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp])*_mag_x[_qp] + 2*_alpha[_qp]*(_mag_x_grad[_qp](0)*_j_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp])*_mag_z[_qp])*_phi[_j][_qp])*_xi[_qp]))/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*_e[_qp]*_g0[_qp]*_Ms[_qp]*(1 + Utility::pow<2>(_xi[_qp])));
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
      return (_muB[_qp]*_test[_i][_qp]*(_grad_phi[_j][_qp](0)*_j_x[_qp]*_mag_x[_qp]*_mag_y[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp]*_mag_x[_qp]*_mag_y[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp]*_mag_x[_qp]*_mag_y[_qp] - _alpha[_qp]*_grad_phi[_j][_qp](0)*_j_x[_qp]*_mag_z[_qp] - _alpha[_qp]*_grad_phi[_j][_qp](1)*_j_y[_qp]*_mag_z[_qp] - _alpha[_qp]*_grad_phi[_j][_qp](2)*_j_z[_qp]*_mag_z[_qp] + _alpha[_qp]*_mag_z_grad[_qp](0)*_j_x[_qp]*_phi[_j][_qp] + _alpha[_qp]*_mag_z_grad[_qp](1)*_j_y[_qp]*_phi[_j][_qp] + _alpha[_qp]*_mag_z_grad[_qp](2)*_j_z[_qp]*_phi[_j][_qp] - 2*_mag_y_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp]*_phi[_j][_qp] - 
       2*_mag_y_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp]*_phi[_j][_qp] - 2*_mag_y_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp]*_phi[_j][_qp] + _mag_x_grad[_qp](0)*_j_x[_qp]*_mag_y[_qp]*_phi[_j][_qp] + _mag_x_grad[_qp](1)*_j_y[_qp]*_mag_y[_qp]*_phi[_j][_qp] + _mag_x_grad[_qp](2)*_j_z[_qp]*_mag_y[_qp]*_phi[_j][_qp] + 
       ((_grad_phi[_j][_qp](0)*_j_x[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp])*(_alpha[_qp]*_mag_x[_qp]*_mag_y[_qp] + _mag_z[_qp]) - (_mag_z_grad[_qp](0)*_j_x[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp] + 2*_alpha[_qp]*(_mag_y_grad[_qp](0)*_j_x[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp])*_mag_x[_qp] - _alpha[_qp]*(_mag_x_grad[_qp](0)*_j_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp])*_mag_y[_qp])*_phi[_j][_qp])*_xi[_qp]))/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*_e[_qp]*_g0[_qp]*_Ms[_qp]*(1 + Utility::pow<2>(_xi[_qp])));
    }
    else if (jvar == _mag_z_var)
    {
      return (_muB[_qp]*_test[_i][_qp]*(_alpha[_qp]*_grad_phi[_j][_qp](0)*_j_x[_qp]*_mag_x[_qp] + _alpha[_qp]*_grad_phi[_j][_qp](1)*_j_y[_qp]*_mag_x[_qp] + _alpha[_qp]*_grad_phi[_j][_qp](2)*_j_z[_qp]*_mag_x[_qp] + _grad_phi[_j][_qp](0)*_j_x[_qp]*_mag_y[_qp]*_mag_z[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp]*_mag_y[_qp]*_mag_z[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp]*_mag_y[_qp]*_mag_z[_qp] - _alpha[_qp]*_mag_x_grad[_qp](0)*_j_x[_qp]*_phi[_j][_qp] - _alpha[_qp]*_mag_x_grad[_qp](1)*_j_y[_qp]*_phi[_j][_qp] - _alpha[_qp]*_mag_x_grad[_qp](2)*_j_z[_qp]*_phi[_j][_qp] + _mag_z_grad[_qp](0)*_j_x[_qp]*_mag_y[_qp]*_phi[_j][_qp] + 
       _mag_z_grad[_qp](1)*_j_y[_qp]*_mag_y[_qp]*_phi[_j][_qp] + _mag_z_grad[_qp](2)*_j_z[_qp]*_mag_y[_qp]*_phi[_j][_qp] - 2*_mag_y_grad[_qp](0)*_j_x[_qp]*_mag_z[_qp]*_phi[_j][_qp] - 2*_mag_y_grad[_qp](1)*_j_y[_qp]*_mag_z[_qp]*_phi[_j][_qp] - 2*_mag_y_grad[_qp](2)*_j_z[_qp]*_mag_z[_qp]*_phi[_j][_qp] + 
       (-((_grad_phi[_j][_qp](0)*_j_x[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp])*(_mag_x[_qp] - _alpha[_qp]*_mag_y[_qp]*_mag_z[_qp])) + (_mag_x_grad[_qp](0)*_j_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp] + _alpha[_qp]*(_mag_z_grad[_qp](0)*_j_x[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp])*_mag_y[_qp] - 2*_alpha[_qp]*(_mag_y_grad[_qp](0)*_j_x[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp])*_mag_z[_qp])*_phi[_j][_qp])*_xi[_qp]))/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*_e[_qp]*_g0[_qp]*_Ms[_qp]*(1 + Utility::pow<2>(_xi[_qp])));
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
      return (_muB[_qp]*_test[_i][_qp]*(_alpha[_qp]*_grad_phi[_j][_qp](0)*_j_x[_qp]*_mag_y[_qp] + _alpha[_qp]*_grad_phi[_j][_qp](1)*_j_y[_qp]*_mag_y[_qp] + _alpha[_qp]*_grad_phi[_j][_qp](2)*_j_z[_qp]*_mag_y[_qp] + _grad_phi[_j][_qp](0)*_j_x[_qp]*_mag_x[_qp]*_mag_z[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp]*_mag_x[_qp]*_mag_z[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp]*_mag_x[_qp]*_mag_z[_qp] - _alpha[_qp]*_mag_y_grad[_qp](0)*_j_x[_qp]*_phi[_j][_qp] - _alpha[_qp]*_mag_y_grad[_qp](1)*_j_y[_qp]*_phi[_j][_qp] - _alpha[_qp]*_mag_y_grad[_qp](2)*_j_z[_qp]*_phi[_j][_qp] - 2*_mag_z_grad[_qp](0)*_j_x[_qp]*_mag_x[_qp]*_phi[_j][_qp] - 
       2*_mag_z_grad[_qp](1)*_j_y[_qp]*_mag_x[_qp]*_phi[_j][_qp] - 2*_mag_z_grad[_qp](2)*_j_z[_qp]*_mag_x[_qp]*_phi[_j][_qp] + _mag_x_grad[_qp](0)*_j_x[_qp]*_mag_z[_qp]*_phi[_j][_qp] + _mag_x_grad[_qp](1)*_j_y[_qp]*_mag_z[_qp]*_phi[_j][_qp] + _mag_x_grad[_qp](2)*_j_z[_qp]*_mag_z[_qp]*_phi[_j][_qp] + 
       (-((_grad_phi[_j][_qp](0)*_j_x[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp])*(_mag_y[_qp] - _alpha[_qp]*_mag_x[_qp]*_mag_z[_qp])) + (_mag_y_grad[_qp](0)*_j_x[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp] - 2*_alpha[_qp]*(_mag_z_grad[_qp](0)*_j_x[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp])*_mag_x[_qp] + _alpha[_qp]*(_mag_x_grad[_qp](0)*_j_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp])*_mag_z[_qp])*_phi[_j][_qp])*_xi[_qp]))/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*_e[_qp]*_g0[_qp]*_Ms[_qp]*(1 + Utility::pow<2>(_xi[_qp])));
    }
    else if (jvar == _mag_y_var)
    {
      return (_muB[_qp]*_test[_i][_qp]*(-(_alpha[_qp]*_grad_phi[_j][_qp](0)*_j_x[_qp]*_mag_x[_qp]) - _alpha[_qp]*_grad_phi[_j][_qp](1)*_j_y[_qp]*_mag_x[_qp] - _alpha[_qp]*_grad_phi[_j][_qp](2)*_j_z[_qp]*_mag_x[_qp] + _grad_phi[_j][_qp](0)*_j_x[_qp]*_mag_y[_qp]*_mag_z[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp]*_mag_y[_qp]*_mag_z[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp]*_mag_y[_qp]*_mag_z[_qp] + _alpha[_qp]*_mag_x_grad[_qp](0)*_j_x[_qp]*_phi[_j][_qp] + _alpha[_qp]*_mag_x_grad[_qp](1)*_j_y[_qp]*_phi[_j][_qp] + _alpha[_qp]*_mag_x_grad[_qp](2)*_j_z[_qp]*_phi[_j][_qp] - 2*_mag_z_grad[_qp](0)*_j_x[_qp]*_mag_y[_qp]*_phi[_j][_qp] - 
       2*_mag_z_grad[_qp](1)*_j_y[_qp]*_mag_y[_qp]*_phi[_j][_qp] - 2*_mag_z_grad[_qp](2)*_j_z[_qp]*_mag_y[_qp]*_phi[_j][_qp] + _mag_y_grad[_qp](0)*_j_x[_qp]*_mag_z[_qp]*_phi[_j][_qp] + _mag_y_grad[_qp](1)*_j_y[_qp]*_mag_z[_qp]*_phi[_j][_qp] + _mag_y_grad[_qp](2)*_j_z[_qp]*_mag_z[_qp]*_phi[_j][_qp] + 
       ((_grad_phi[_j][_qp](0)*_j_x[_qp] + _grad_phi[_j][_qp](1)*_j_y[_qp] + _grad_phi[_j][_qp](2)*_j_z[_qp])*(_mag_x[_qp] + _alpha[_qp]*_mag_y[_qp]*_mag_z[_qp]) - (_mag_x_grad[_qp](0)*_j_x[_qp] + _mag_x_grad[_qp](1)*_j_y[_qp] + _mag_x_grad[_qp](2)*_j_z[_qp] + 2*_alpha[_qp]*(_mag_z_grad[_qp](0)*_j_x[_qp] + _mag_z_grad[_qp](1)*_j_y[_qp] + _mag_z_grad[_qp](2)*_j_z[_qp])*_mag_y[_qp] - _alpha[_qp]*(_mag_y_grad[_qp](0)*_j_x[_qp] + _mag_y_grad[_qp](1)*_j_y[_qp] + _mag_y_grad[_qp](2)*_j_z[_qp])*_mag_z[_qp])*_phi[_j][_qp])*_xi[_qp]))/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*_e[_qp]*_g0[_qp]*_Ms[_qp]*(1 + Utility::pow<2>(_xi[_qp])));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
