/**
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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "DampingMagneticExchangeDerivative.h"

template<>
InputParameters validParams<DampingMagneticExchangeDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates the residual for exchange term in the damped part of the LLG equation.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the antiferromagnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the antiferromagnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the antiferromagnetic vector");
  params.addRequiredParam<Real>("A", "The exchange coefficient");
  params.addRequiredParam<Real>("M0", "The scaled magnetization vector value");
  params.addRequiredParam<Real>("alphaLL", "The damping coefficient");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

DampingMagneticExchangeDerivative::DampingMagneticExchangeDerivative(const InputParameters & parameters)
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
  _A(getParam<Real>("A")),
  _M0(getParam<Real>("M0")),
  _alphaLL(getParam<Real>("alphaLL")),
  _len_scale(getParam<Real>("len_scale"))
{
}


Real
DampingMagneticExchangeDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    RealVectorValue w(-(_mag_y_grad[_qp](0)*_mag_x[_qp]*_mag_y[_qp]) + _mag_x_grad[_qp](0)*std::pow(_mag_y[_qp],2) - _mag_z_grad[_qp](0)*_mag_x[_qp]*_mag_z[_qp] + _mag_x_grad[_qp](0)*std::pow(_mag_z[_qp],2), -(_mag_y_grad[_qp](1)*_mag_x[_qp]*_mag_y[_qp]) + _mag_x_grad[_qp](1)*std::pow(_mag_y[_qp],2) - _mag_z_grad[_qp](1)*_mag_x[_qp]*_mag_z[_qp] + _mag_x_grad[_qp](1)*std::pow(_mag_z[_qp],2), -(_mag_y_grad[_qp](2)*_mag_x[_qp]*_mag_y[_qp]) + _mag_x_grad[_qp](2)*std::pow(_mag_y[_qp],2) - _mag_z_grad[_qp](2)*_mag_x[_qp]*_mag_z[_qp] + _mag_x_grad[_qp](2)*std::pow(_mag_z[_qp],2));
    return _alphaLL * 2.0 * _A * (1.0 / (8 * _M0 *_M0 * _M0)) * ( _grad_test[_i][_qp] * w );
  }
  else if (_component == 1)
  {
    RealVectorValue w(_mag_y_grad[_qp](0)*std::pow(_mag_x[_qp],2) - _mag_x_grad[_qp](0)*_mag_x[_qp]*_mag_y[_qp] - _mag_z_grad[_qp](0)*_mag_y[_qp]*_mag_z[_qp] + _mag_y_grad[_qp](0)*std::pow(_mag_z[_qp],2), _mag_y_grad[_qp](1)*std::pow(_mag_x[_qp],2) - _mag_x_grad[_qp](1)*_mag_x[_qp]*_mag_y[_qp] - _mag_z_grad[_qp](1)*_mag_y[_qp]*_mag_z[_qp] + _mag_y_grad[_qp](1)*std::pow(_mag_z[_qp],2), _mag_y_grad[_qp](2)*std::pow(_mag_x[_qp],2) - _mag_x_grad[_qp](2)*_mag_x[_qp]*_mag_y[_qp] - _mag_z_grad[_qp](2)*_mag_y[_qp]*_mag_z[_qp] + _mag_y_grad[_qp](2)*std::pow(_mag_z[_qp],2));
    return _alphaLL * 2.0 * _A * (1.0 / (8 * _M0 *_M0 * _M0)) * ( _grad_test[_i][_qp] * w );
  }
  else if (_component == 2)
  {
    RealVectorValue w(_mag_z_grad[_qp](0)*std::pow(_mag_x[_qp],2) + _mag_z_grad[_qp](0)*std::pow(_mag_y[_qp],2) - _mag_x_grad[_qp](0)*_mag_x[_qp]*_mag_z[_qp] - _mag_y_grad[_qp](0)*_mag_y[_qp]*_mag_z[_qp], _mag_z_grad[_qp](1)*std::pow(_mag_x[_qp],2) + _mag_z_grad[_qp](1)*std::pow(_mag_y[_qp],2) - _mag_x_grad[_qp](1)*_mag_x[_qp]*_mag_z[_qp] - _mag_y_grad[_qp](1)*_mag_y[_qp]*_mag_z[_qp], _mag_z_grad[_qp](2)*std::pow(_mag_x[_qp],2) + _mag_z_grad[_qp](2)*std::pow(_mag_y[_qp],2) - _mag_x_grad[_qp](2)*_mag_x[_qp]*_mag_z[_qp] - _mag_y_grad[_qp](2)*_mag_y[_qp]*_mag_z[_qp]);
    return _alphaLL * 2.0 * _A * (1.0 / (8 * _M0 *_M0 * _M0)) * ( _grad_test[_i][_qp] * w );
  }
  else
    return 0.0;
}

Real
DampingMagneticExchangeDerivative::computeQpJacobian()
{
  if (_component == 0)
  {
    RealVectorValue w(-(_mag_y_grad[_qp](0)*_mag_y[_qp]) - _mag_z_grad[_qp](0)*_mag_z[_qp], -(_mag_y_grad[_qp](1)*_mag_y[_qp]) - _mag_z_grad[_qp](1)*_mag_z[_qp], -(_mag_y_grad[_qp](2)*_mag_y[_qp]) - _mag_z_grad[_qp](2)*_mag_z[_qp]);
    return _alphaLL * 2.0 * _A * (1.0 / (8 * _M0 *_M0 * _M0)) * _phi[_j][_qp] * ( _grad_test[_i][_qp] * w );
  }
  else if (_component == 1)
  {
    RealVectorValue w(-(_mag_x_grad[_qp](0)*_mag_x[_qp]) - _mag_z_grad[_qp](0)*_mag_z[_qp], -(_mag_x_grad[_qp](1)*_mag_x[_qp]) - _mag_z_grad[_qp](1)*_mag_z[_qp], -(_mag_x_grad[_qp](2)*_mag_x[_qp]) - _mag_z_grad[_qp](2)*_mag_z[_qp]);
    return _alphaLL * 2.0 * _A * (1.0 / (8 * _M0 *_M0 * _M0)) * _phi[_j][_qp] * ( _grad_test[_i][_qp] * w );
  }
  else if (_component == 2)
  {
    RealVectorValue w(-(_mag_x_grad[_qp](0)*_mag_x[_qp]) - _mag_y_grad[_qp](0)*_mag_y[_qp], -(_mag_x_grad[_qp](1)*_mag_x[_qp]) - _mag_y_grad[_qp](1)*_mag_y[_qp], -(_mag_x_grad[_qp](2)*_mag_x[_qp]) - _mag_y_grad[_qp](2)*_mag_y[_qp]);
    return _alphaLL * 2.0 * _A * (1.0 / (8 * _M0 *_M0 * _M0)) * _phi[_j][_qp] * ( _grad_test[_i][_qp] * w );
  }
  else
    return 0.0;
}

Real
DampingMagneticExchangeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      RealVectorValue w(-(_grad_phi[_j][_qp](0)*_mag_x[_qp]*_mag_y[_qp]) + _phi[_j][_qp]*(-(_mag_y_grad[_qp](0)*_mag_x[_qp]) + 2*_mag_x_grad[_qp](0)*_mag_y[_qp]), -(_grad_phi[_j][_qp](1)*_mag_x[_qp]*_mag_y[_qp]) + _phi[_j][_qp]*(-(_mag_y_grad[_qp](1)*_mag_x[_qp]) + 2*_mag_x_grad[_qp](1)*_mag_y[_qp]), -(_grad_phi[_j][_qp](2)*_mag_x[_qp]*_mag_y[_qp]) + _phi[_j][_qp]*(-(_mag_y_grad[_qp](2)*_mag_x[_qp]) + 2*_mag_x_grad[_qp](2)*_mag_y[_qp]));
      return _alphaLL * 2.0 * _A * (1.0 / (8 * _M0 *_M0 * _M0)) * ( _grad_test[_i][_qp] * w );
    }
    else if (jvar == _mag_z_var)
      {
        RealVectorValue w(-(_grad_phi[_j][_qp](0)*_mag_x[_qp]*_mag_z[_qp]) + _phi[_j][_qp]*(-(_mag_z_grad[_qp](0)*_mag_x[_qp]) + 2*_mag_x_grad[_qp](0)*_mag_z[_qp]), -(_grad_phi[_j][_qp](1)*_mag_x[_qp]*_mag_z[_qp]) + _phi[_j][_qp]*(-(_mag_z_grad[_qp](1)*_mag_x[_qp]) + 2*_mag_x_grad[_qp](1)*_mag_z[_qp]), -(_grad_phi[_j][_qp](2)*_mag_x[_qp]*_mag_z[_qp]) + _phi[_j][_qp]*(-(_mag_z_grad[_qp](2)*_mag_x[_qp]) + 2*_mag_x_grad[_qp](2)*_mag_z[_qp]));
        return _alphaLL * 2.0 * _A * (1.0 / (8 * _M0 *_M0 * _M0)) * ( _grad_test[_i][_qp] * w );
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
        RealVectorValue w(-(_grad_phi[_j][_qp](0)*_mag_x[_qp]*_mag_y[_qp]) + _phi[_j][_qp]*(2*_mag_y_grad[_qp](0)*_mag_x[_qp] - _mag_x_grad[_qp](0)*_mag_y[_qp]), -(_grad_phi[_j][_qp](1)*_mag_x[_qp]*_mag_y[_qp]) + _phi[_j][_qp]*(2*_mag_y_grad[_qp](1)*_mag_x[_qp] - _mag_x_grad[_qp](1)*_mag_y[_qp]), -(_grad_phi[_j][_qp](2)*_mag_x[_qp]*_mag_y[_qp]) + _phi[_j][_qp]*(2*_mag_y_grad[_qp](2)*_mag_x[_qp] - _mag_x_grad[_qp](2)*_mag_y[_qp]));
        return _alphaLL * 2.0 * _A * (1.0 / (8 * _M0 *_M0 * _M0)) * ( _grad_test[_i][_qp] * w );
      }
    else if (jvar == _mag_z_var)
      {
        RealVectorValue w(-(_grad_phi[_j][_qp](0)*_mag_y[_qp]*_mag_z[_qp]) + _phi[_j][_qp]*(-(_mag_z_grad[_qp](0)*_mag_y[_qp]) + 2*_mag_y_grad[_qp](0)*_mag_z[_qp]), -(_grad_phi[_j][_qp](1)*_mag_y[_qp]*_mag_z[_qp]) + _phi[_j][_qp]*(-(_mag_z_grad[_qp](1)*_mag_y[_qp]) + 2*_mag_y_grad[_qp](1)*_mag_z[_qp]), -(_grad_phi[_j][_qp](2)*_mag_y[_qp]*_mag_z[_qp]) + _phi[_j][_qp]*(-(_mag_z_grad[_qp](2)*_mag_y[_qp]) + 2*_mag_y_grad[_qp](2)*_mag_z[_qp]));
        return _alphaLL * 2.0 * _A * (1.0 / (8 * _M0 *_M0 * _M0)) * ( _grad_test[_i][_qp] * w );
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
        RealVectorValue w(-(_grad_phi[_j][_qp](0)*_mag_x[_qp]*_mag_z[_qp]) + _phi[_j][_qp]*(2*_mag_z_grad[_qp](0)*_mag_x[_qp] - _mag_x_grad[_qp](0)*_mag_z[_qp]), -(_grad_phi[_j][_qp](1)*_mag_x[_qp]*_mag_z[_qp]) + _phi[_j][_qp]*(2*_mag_z_grad[_qp](1)*_mag_x[_qp] - _mag_x_grad[_qp](1)*_mag_z[_qp]), -(_grad_phi[_j][_qp](2)*_mag_x[_qp]*_mag_z[_qp]) + _phi[_j][_qp]*(2*_mag_z_grad[_qp](2)*_mag_x[_qp] - _mag_x_grad[_qp](2)*_mag_z[_qp]));
        return _alphaLL * 2.0 * _A * (1.0 / (8 * _M0 *_M0 * _M0)) * ( _grad_test[_i][_qp] * w );
      }
    else if (jvar == _mag_y_var)
      {
        RealVectorValue w(-(_grad_phi[_j][_qp](0)*_mag_y[_qp]*_mag_z[_qp]) + _phi[_j][_qp]*(2*_mag_z_grad[_qp](0)*_mag_y[_qp] - _mag_y_grad[_qp](0)*_mag_z[_qp]), -(_grad_phi[_j][_qp](1)*_mag_y[_qp]*_mag_z[_qp]) + _phi[_j][_qp]*(2*_mag_z_grad[_qp](1)*_mag_y[_qp] - _mag_y_grad[_qp](1)*_mag_z[_qp]), -(_grad_phi[_j][_qp](2)*_mag_y[_qp]*_mag_z[_qp]) + _phi[_j][_qp]*(2*_mag_z_grad[_qp](2)*_mag_y[_qp] - _mag_y_grad[_qp](2)*_mag_z[_qp]));
        return _alphaLL * 2.0 * _A * (1.0 / (8 * _M0 *_M0 * _M0)) * ( _grad_test[_i][_qp] * w );
      }
    else
      {
        return 0.0;
      }
  }
  else
    return 0.0;
}
