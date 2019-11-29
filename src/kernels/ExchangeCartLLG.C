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

#include "ExchangeCartLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", ExchangeCartLLG);

template<>
InputParameters validParams<ExchangeCartLLG>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic anisotropy energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  params.addRequiredParam<Real>("alpha", "the damping coefficient in the LLG equation");
  params.addRequiredParam<Real>("g0", "g0");
  params.addRequiredParam<Real>("Ae", "Ae");
  params.addRequiredParam<Real>("Ms", "Ms");
  return params;
}

ExchangeCartLLG::ExchangeCartLLG(const InputParameters & parameters)
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
  _alpha(getParam<Real>("alpha")),
  _g0(getParam<Real>("g0")),
  _Ae(getParam<Real>("Ae")),
  _Ms(getParam<Real>("Ms"))
{
}

Real
ExchangeCartLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return (-2*_Ae*_g0*(_alpha*(_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_x_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_x_grad[_qp](2)*_grad_test[_i][_qp](2))*(-1 + Utility::pow<2>(_mag_x[_qp])) + (_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_y_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_y_grad[_qp](2)*_grad_test[_i][_qp](2))*(_mag_x[_qp]*_mag_y[_qp] - _mag_z[_qp]) + 
       (_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_z_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_z_grad[_qp](2)*_grad_test[_i][_qp](2))*(_mag_y[_qp] + _alpha*_mag_x[_qp]*_mag_z[_qp]) + 2*_alpha*(Utility::pow<2>(_mag_x_grad[_qp](0)) + Utility::pow<2>(_mag_x_grad[_qp](1)) + Utility::pow<2>(_mag_x_grad[_qp](2)))*_mag_x[_qp]*_test[_i][_qp] + 
       (_mag_y_grad[_qp](0)*(-_mag_z_grad[_qp](0) + _mag_y_grad[_qp](0)*_mag_x[_qp] + _mag_x_grad[_qp](0)*_mag_y[_qp]) + _mag_y_grad[_qp](1)*(-_mag_z_grad[_qp](1) + _mag_y_grad[_qp](1)*_mag_x[_qp] + _mag_x_grad[_qp](1)*_mag_y[_qp]) + _mag_y_grad[_qp](2)*(-_mag_z_grad[_qp](2) + _mag_y_grad[_qp](2)*_mag_x[_qp] + _mag_x_grad[_qp](2)*_mag_y[_qp]))*_test[_i][_qp] + 
       (_mag_y_grad[_qp](0)*_mag_z_grad[_qp](0) + _mag_y_grad[_qp](1)*_mag_z_grad[_qp](1) + _mag_y_grad[_qp](2)*_mag_z_grad[_qp](2) + _alpha*(Utility::pow<2>(_mag_z_grad[_qp](0)) + Utility::pow<2>(_mag_z_grad[_qp](1)) + Utility::pow<2>(_mag_z_grad[_qp](2)))*_mag_x[_qp] + _alpha*(_mag_x_grad[_qp](0)*_mag_z_grad[_qp](0) + _mag_x_grad[_qp](1)*_mag_z_grad[_qp](1) + _mag_x_grad[_qp](2)*_mag_z_grad[_qp](2))*_mag_z[_qp])*_test[_i][_qp]))/((1 + Utility::pow<2>(_alpha))*_Ms);
  }
  else if (_component == 1)
  {
    return (-2*_Ae*_g0*(-((_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_z_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_z_grad[_qp](2)*_grad_test[_i][_qp](2))*_mag_x[_qp]) + (_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_x_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_x_grad[_qp](2)*_grad_test[_i][_qp](2))*_mag_z[_qp] + 
       _alpha*(-(_mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)) - _mag_y_grad[_qp](2)*_grad_test[_i][_qp](2) + _mag_y[_qp]*(_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp] + _mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp] + _mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp] + _mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp] + _mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp] + _mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_z[_qp] + _mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_z[_qp] + _mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_z[_qp]) + 
          2*Utility::pow<2>(_mag_y_grad[_qp](0))*_mag_y[_qp]*_test[_i][_qp] + (_mag_x_grad[_qp](1)*_mag_y_grad[_qp](1)*_mag_x[_qp] + _mag_x_grad[_qp](2)*_mag_y_grad[_qp](2)*_mag_x[_qp] + Utility::pow<2>(_mag_x_grad[_qp](1))*_mag_y[_qp] + Utility::pow<2>(_mag_x_grad[_qp](2))*_mag_y[_qp] + 
             (Utility::pow<2>(_mag_x_grad[_qp](0)) + 2*Utility::pow<2>(_mag_y_grad[_qp](1)) + 2*Utility::pow<2>(_mag_y_grad[_qp](2)) + Utility::pow<2>(_mag_z_grad[_qp](0)) + Utility::pow<2>(_mag_z_grad[_qp](1)) + Utility::pow<2>(_mag_z_grad[_qp](2)))*_mag_y[_qp] + (_mag_y_grad[_qp](1)*_mag_z_grad[_qp](1) + _mag_y_grad[_qp](2)*_mag_z_grad[_qp](2))*_mag_z[_qp])*_test[_i][_qp] + 
          _mag_y_grad[_qp](0)*(_grad_test[_i][_qp](0)*(-1 + Utility::pow<2>(_mag_y[_qp])) + (_mag_x_grad[_qp](0)*_mag_x[_qp] + _mag_z_grad[_qp](0)*_mag_z[_qp])*_test[_i][_qp]))))/((1 + Utility::pow<2>(_alpha))*_Ms);
  }
  else if (_component == 2)
  {
    return (-2*_Ae*_g0*((_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_y_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_y_grad[_qp](2)*_grad_test[_i][_qp](2))*_mag_x[_qp] - (_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_x_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_x_grad[_qp](2)*_grad_test[_i][_qp](2))*_mag_y[_qp] + 
       _alpha*(-(_mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)) - _mag_z_grad[_qp](2)*_grad_test[_i][_qp](2) + _mag_z[_qp]*(_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp] + _mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp] + _mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp] + _mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp] + _mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp] + _mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp] + _mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_z[_qp] + _mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_z[_qp]) + 
          2*Utility::pow<2>(_mag_z_grad[_qp](0))*_mag_z[_qp]*_test[_i][_qp] + (_mag_x_grad[_qp](1)*_mag_z_grad[_qp](1)*_mag_x[_qp] + _mag_x_grad[_qp](2)*_mag_z_grad[_qp](2)*_mag_x[_qp] + _mag_y_grad[_qp](1)*_mag_z_grad[_qp](1)*_mag_y[_qp] + _mag_y_grad[_qp](2)*_mag_z_grad[_qp](2)*_mag_y[_qp] + 
             (Utility::pow<2>(_mag_x_grad[_qp](0)) + Utility::pow<2>(_mag_x_grad[_qp](1)) + Utility::pow<2>(_mag_x_grad[_qp](2)) + Utility::pow<2>(_mag_y_grad[_qp](0)) + Utility::pow<2>(_mag_y_grad[_qp](1)) + Utility::pow<2>(_mag_y_grad[_qp](2)) + 2*(Utility::pow<2>(_mag_z_grad[_qp](1)) + Utility::pow<2>(_mag_z_grad[_qp](2))))*_mag_z[_qp])*_test[_i][_qp] + 
          _mag_z_grad[_qp](0)*(_grad_test[_i][_qp](0)*(-1 + Utility::pow<2>(_mag_z[_qp])) + (_mag_x_grad[_qp](0)*_mag_x[_qp] + _mag_y_grad[_qp](0)*_mag_y[_qp])*_test[_i][_qp]))))/((1 + Utility::pow<2>(_alpha))*_Ms);
  }
  else
    return 0.0;
}

Real
ExchangeCartLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return (-2*_Ae*_g0*(_alpha*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*(-1 + Utility::pow<2>(_mag_x[_qp])) + 2*_alpha*(_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_x_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_x_grad[_qp](2)*_grad_test[_i][_qp](2))*_mag_x[_qp]*_phi[_j][_qp] + (_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_y_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_y_grad[_qp](2)*_grad_test[_i][_qp](2))*_mag_y[_qp]*_phi[_j][_qp] + 
       _alpha*(_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_z_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_z_grad[_qp](2)*_grad_test[_i][_qp](2))*_mag_z[_qp]*_phi[_j][_qp] + 4*_alpha*(_mag_x_grad[_qp](0)*_grad_phi[_j][_qp](0) + _mag_x_grad[_qp](1)*_grad_phi[_j][_qp](1) + _mag_x_grad[_qp](2)*_grad_phi[_j][_qp](2))*_mag_x[_qp]*_test[_i][_qp] + 2*_alpha*(Utility::pow<2>(_mag_x_grad[_qp](0)) + Utility::pow<2>(_mag_x_grad[_qp](1)) + Utility::pow<2>(_mag_x_grad[_qp](2)))*_phi[_j][_qp]*_test[_i][_qp] + 
       ((_mag_y_grad[_qp](0)*_grad_phi[_j][_qp](0) + _mag_y_grad[_qp](1)*_grad_phi[_j][_qp](1) + _mag_y_grad[_qp](2)*_grad_phi[_j][_qp](2))*_mag_y[_qp] + (Utility::pow<2>(_mag_y_grad[_qp](0)) + Utility::pow<2>(_mag_y_grad[_qp](1)) + Utility::pow<2>(_mag_y_grad[_qp](2)))*_phi[_j][_qp])*_test[_i][_qp] + 
       (_alpha*(_mag_z_grad[_qp](0)*_grad_phi[_j][_qp](0) + _mag_z_grad[_qp](1)*_grad_phi[_j][_qp](1) + _mag_z_grad[_qp](2)*_grad_phi[_j][_qp](2))*_mag_z[_qp] + _alpha*(Utility::pow<2>(_mag_z_grad[_qp](0)) + Utility::pow<2>(_mag_z_grad[_qp](1)) + Utility::pow<2>(_mag_z_grad[_qp](2)))*_phi[_j][_qp])*_test[_i][_qp]))/((1 + Utility::pow<2>(_alpha))*_Ms);
  }
  else if (_component == 1)
  {
    return (-2*_Ae*_alpha*_g0*(-(_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_y[_qp]) + _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*(-1 + Utility::pow<2>(_mag_y[_qp])) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*(-1 + Utility::pow<2>(_mag_y[_qp])) + _mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp]*_phi[_j][_qp] + _mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp]*_phi[_j][_qp] + _mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp]*_phi[_j][_qp] + 
       2*_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp]*_phi[_j][_qp] + 2*_mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp]*_phi[_j][_qp] + 2*_mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp]*_phi[_j][_qp] + _mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_z[_qp]*_phi[_j][_qp] + _mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_z[_qp]*_phi[_j][_qp] + _mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_z[_qp]*_phi[_j][_qp] + _grad_phi[_j][_qp](0)*(_mag_x_grad[_qp](0)*_mag_x[_qp] + 4*_mag_y_grad[_qp](0)*_mag_y[_qp] + _mag_z_grad[_qp](0)*_mag_z[_qp])*_test[_i][_qp] + 
       _grad_phi[_j][_qp](1)*(_mag_x_grad[_qp](1)*_mag_x[_qp] + 4*_mag_y_grad[_qp](1)*_mag_y[_qp] + _mag_z_grad[_qp](1)*_mag_z[_qp])*_test[_i][_qp] + (_grad_phi[_j][_qp](2)*(_mag_x_grad[_qp](2)*_mag_x[_qp] + 4*_mag_y_grad[_qp](2)*_mag_y[_qp] + _mag_z_grad[_qp](2)*_mag_z[_qp]) + 
          (Utility::pow<2>(_mag_x_grad[_qp](0)) + Utility::pow<2>(_mag_x_grad[_qp](1)) + Utility::pow<2>(_mag_x_grad[_qp](2)) + 2*Utility::pow<2>(_mag_y_grad[_qp](0)) + 2*Utility::pow<2>(_mag_y_grad[_qp](1)) + 2*Utility::pow<2>(_mag_y_grad[_qp](2)) + Utility::pow<2>(_mag_z_grad[_qp](0)) + Utility::pow<2>(_mag_z_grad[_qp](1)) + Utility::pow<2>(_mag_z_grad[_qp](2)))*_phi[_j][_qp])*_test[_i][_qp]))/((1 + Utility::pow<2>(_alpha))*_Ms);
  }
  else if (_component == 2)
  {
    return (-2*_Ae*_alpha*_g0*(-(_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_z[_qp]) + _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*(-1 + Utility::pow<2>(_mag_z[_qp])) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*(-1 + Utility::pow<2>(_mag_z[_qp])) + _mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp]*_phi[_j][_qp] + _mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp]*_phi[_j][_qp] + _mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp]*_phi[_j][_qp] + 
       _mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp]*_phi[_j][_qp] + _mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp]*_phi[_j][_qp] + _mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp]*_phi[_j][_qp] + 2*_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_z[_qp]*_phi[_j][_qp] + 2*_mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_z[_qp]*_phi[_j][_qp] + 2*_mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_z[_qp]*_phi[_j][_qp] + _grad_phi[_j][_qp](0)*(_mag_x_grad[_qp](0)*_mag_x[_qp] + _mag_y_grad[_qp](0)*_mag_y[_qp] + 4*_mag_z_grad[_qp](0)*_mag_z[_qp])*_test[_i][_qp] + 
       _grad_phi[_j][_qp](1)*(_mag_x_grad[_qp](1)*_mag_x[_qp] + _mag_y_grad[_qp](1)*_mag_y[_qp] + 4*_mag_z_grad[_qp](1)*_mag_z[_qp])*_test[_i][_qp] + (_grad_phi[_j][_qp](2)*(_mag_x_grad[_qp](2)*_mag_x[_qp] + _mag_y_grad[_qp](2)*_mag_y[_qp] + 4*_mag_z_grad[_qp](2)*_mag_z[_qp]) + 
          (Utility::pow<2>(_mag_x_grad[_qp](0)) + Utility::pow<2>(_mag_x_grad[_qp](1)) + Utility::pow<2>(_mag_x_grad[_qp](2)) + Utility::pow<2>(_mag_y_grad[_qp](0)) + Utility::pow<2>(_mag_y_grad[_qp](1)) + Utility::pow<2>(_mag_y_grad[_qp](2)) + 2*(Utility::pow<2>(_mag_z_grad[_qp](0)) + Utility::pow<2>(_mag_z_grad[_qp](1)) + Utility::pow<2>(_mag_z_grad[_qp](2))))*_phi[_j][_qp])*_test[_i][_qp]))/((1 + Utility::pow<2>(_alpha))*_Ms);
  }
  else
    return 0.0;
}

Real
ExchangeCartLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return (2*_Ae*_g0*(-((_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*(_mag_x[_qp]*_mag_y[_qp] - _mag_z[_qp])) - (_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_z_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_z_grad[_qp](2)*_grad_test[_i][_qp](2) + _mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp] + _mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp] + _mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp])*_phi[_j][_qp] - 
       (2*_mag_y_grad[_qp](0)*_grad_phi[_j][_qp](0)*_mag_x[_qp] + 2*_mag_y_grad[_qp](1)*_grad_phi[_j][_qp](1)*_mag_x[_qp] + 2*_mag_y_grad[_qp](2)*_grad_phi[_j][_qp](2)*_mag_x[_qp] + _mag_x_grad[_qp](0)*_grad_phi[_j][_qp](0)*_mag_y[_qp] + _mag_x_grad[_qp](1)*_grad_phi[_j][_qp](1)*_mag_y[_qp] + _mag_x_grad[_qp](2)*_grad_phi[_j][_qp](2)*_mag_y[_qp] + _mag_x_grad[_qp](0)*_mag_y_grad[_qp](0)*_phi[_j][_qp] + _mag_x_grad[_qp](1)*_mag_y_grad[_qp](1)*_phi[_j][_qp] + _mag_x_grad[_qp](2)*_mag_y_grad[_qp](2)*_phi[_j][_qp])*_test[_i][_qp]))/((1 + Utility::pow<2>(_alpha))*_Ms);
    }
    else if (jvar == _mag_z_var)
    {
      return (2*_Ae*_g0*(-((_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*(_mag_y[_qp] + _alpha*_mag_x[_qp]*_mag_z[_qp])) + (_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_y_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_y_grad[_qp](2)*_grad_test[_i][_qp](2) - _alpha*(_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_z_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_z_grad[_qp](2)*_grad_test[_i][_qp](2))*_mag_x[_qp])*_phi[_j][_qp] - 
       _alpha*(2*_mag_z_grad[_qp](0)*_grad_phi[_j][_qp](0)*_mag_x[_qp] + 2*_mag_z_grad[_qp](1)*_grad_phi[_j][_qp](1)*_mag_x[_qp] + 2*_mag_z_grad[_qp](2)*_grad_phi[_j][_qp](2)*_mag_x[_qp] + _mag_x_grad[_qp](0)*_grad_phi[_j][_qp](0)*_mag_z[_qp] + _mag_x_grad[_qp](1)*_grad_phi[_j][_qp](1)*_mag_z[_qp] + _mag_x_grad[_qp](2)*_grad_phi[_j][_qp](2)*_mag_z[_qp] + _mag_x_grad[_qp](0)*_mag_z_grad[_qp](0)*_phi[_j][_qp] + _mag_x_grad[_qp](1)*_mag_z_grad[_qp](1)*_phi[_j][_qp] + _mag_x_grad[_qp](2)*_mag_z_grad[_qp](2)*_phi[_j][_qp])*_test[_i][_qp]))/((1 + Utility::pow<2>(_alpha))*_Ms);
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
      return (2*_Ae*_g0*(-((_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*(_alpha*_mag_x[_qp]*_mag_y[_qp] + _mag_z[_qp])) + (_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_z_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_z_grad[_qp](2)*_grad_test[_i][_qp](2) - _alpha*(_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_x_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_x_grad[_qp](2)*_grad_test[_i][_qp](2))*_mag_y[_qp])*_phi[_j][_qp] - 
       _alpha*(_mag_y_grad[_qp](0)*_grad_phi[_j][_qp](0)*_mag_x[_qp] + _mag_y_grad[_qp](1)*_grad_phi[_j][_qp](1)*_mag_x[_qp] + _mag_y_grad[_qp](2)*_grad_phi[_j][_qp](2)*_mag_x[_qp] + 2*_mag_x_grad[_qp](0)*_grad_phi[_j][_qp](0)*_mag_y[_qp] + 2*_mag_x_grad[_qp](1)*_grad_phi[_j][_qp](1)*_mag_y[_qp] + 2*_mag_x_grad[_qp](2)*_grad_phi[_j][_qp](2)*_mag_y[_qp] + _mag_x_grad[_qp](0)*_mag_y_grad[_qp](0)*_phi[_j][_qp] + _mag_x_grad[_qp](1)*_mag_y_grad[_qp](1)*_phi[_j][_qp] + _mag_x_grad[_qp](2)*_mag_y_grad[_qp](2)*_phi[_j][_qp])*_test[_i][_qp]))/((1 + Utility::pow<2>(_alpha))*_Ms);
    }
    else if (jvar == _mag_z_var)
    {
      return (-2*_Ae*_g0*(-((_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*_mag_x[_qp]) + (_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_x_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_x_grad[_qp](2)*_grad_test[_i][_qp](2))*_phi[_j][_qp] + 
       _alpha*_mag_y[_qp]*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*_mag_z[_qp] + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*_mag_z[_qp] + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*_mag_z[_qp] + _mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*_phi[_j][_qp] + _mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_phi[_j][_qp] + _mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_phi[_j][_qp]) + 
       _alpha*(2*_mag_z_grad[_qp](0)*_grad_phi[_j][_qp](0)*_mag_y[_qp] + 2*_mag_z_grad[_qp](1)*_grad_phi[_j][_qp](1)*_mag_y[_qp] + 2*_mag_z_grad[_qp](2)*_grad_phi[_j][_qp](2)*_mag_y[_qp] + _mag_y_grad[_qp](0)*_grad_phi[_j][_qp](0)*_mag_z[_qp] + _mag_y_grad[_qp](1)*_grad_phi[_j][_qp](1)*_mag_z[_qp] + _mag_y_grad[_qp](2)*_grad_phi[_j][_qp](2)*_mag_z[_qp] + _mag_y_grad[_qp](0)*_mag_z_grad[_qp](0)*_phi[_j][_qp] + _mag_y_grad[_qp](1)*_mag_z_grad[_qp](1)*_phi[_j][_qp] + _mag_y_grad[_qp](2)*_mag_z_grad[_qp](2)*_phi[_j][_qp])*_test[_i][_qp]))/((1 + Utility::pow<2>(_alpha))*_Ms);
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
      return (-2*_Ae*_g0*(-((_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*_mag_y[_qp]) + (_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_y_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_y_grad[_qp](2)*_grad_test[_i][_qp](2))*_phi[_j][_qp] + 
       _alpha*_mag_z[_qp]*(_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp] + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp] + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp] + _mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_phi[_j][_qp] + _mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_phi[_j][_qp] + _mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_phi[_j][_qp]) + 
       _alpha*(_mag_z_grad[_qp](0)*_grad_phi[_j][_qp](0)*_mag_x[_qp] + _mag_z_grad[_qp](1)*_grad_phi[_j][_qp](1)*_mag_x[_qp] + _mag_z_grad[_qp](2)*_grad_phi[_j][_qp](2)*_mag_x[_qp] + 2*_mag_x_grad[_qp](0)*_grad_phi[_j][_qp](0)*_mag_z[_qp] + 2*_mag_x_grad[_qp](1)*_grad_phi[_j][_qp](1)*_mag_z[_qp] + 2*_mag_x_grad[_qp](2)*_grad_phi[_j][_qp](2)*_mag_z[_qp] + _mag_x_grad[_qp](0)*_mag_z_grad[_qp](0)*_phi[_j][_qp] + _mag_x_grad[_qp](1)*_mag_z_grad[_qp](1)*_phi[_j][_qp] + _mag_x_grad[_qp](2)*_mag_z_grad[_qp](2)*_phi[_j][_qp])*_test[_i][_qp]))/((1 + Utility::pow<2>(_alpha))*_Ms);
    }
    else if (jvar == _mag_y_var)
    {
      return (2*_Ae*_g0*(-((_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*(_mag_x[_qp] + _alpha*_mag_y[_qp]*_mag_z[_qp])) + (_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_x_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_x_grad[_qp](2)*_grad_test[_i][_qp](2) - _alpha*(_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0) + _mag_y_grad[_qp](1)*_grad_test[_i][_qp](1) + _mag_y_grad[_qp](2)*_grad_test[_i][_qp](2))*_mag_z[_qp])*_phi[_j][_qp] - 
       _alpha*(_mag_z_grad[_qp](0)*_grad_phi[_j][_qp](0)*_mag_y[_qp] + _mag_z_grad[_qp](1)*_grad_phi[_j][_qp](1)*_mag_y[_qp] + _mag_z_grad[_qp](2)*_grad_phi[_j][_qp](2)*_mag_y[_qp] + 2*_mag_y_grad[_qp](0)*_grad_phi[_j][_qp](0)*_mag_z[_qp] + 2*_mag_y_grad[_qp](1)*_grad_phi[_j][_qp](1)*_mag_z[_qp] + 2*_mag_y_grad[_qp](2)*_grad_phi[_j][_qp](2)*_mag_z[_qp] + _mag_y_grad[_qp](0)*_mag_z_grad[_qp](0)*_phi[_j][_qp] + _mag_y_grad[_qp](1)*_mag_z_grad[_qp](1)*_phi[_j][_qp] + _mag_y_grad[_qp](2)*_mag_z_grad[_qp](2)*_phi[_j][_qp])*_test[_i][_qp]))/((1 + Utility::pow<2>(_alpha))*_Ms);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
