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

#include "SpinFlexoelectricPolarDerivative.h"

template<>
InputParameters validParams<SpinFlexoelectricPolarDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization vector");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization vector");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization vector");
  params.addRequiredCoupledVar("antiferromag_L_x", "The x component of the antiferromagnetic vector");
  params.addCoupledVar("antiferromag_L_y", 0.0, "The y component of the antiferromagnetic vector");
  params.addCoupledVar("antiferromag_L_z", 0.0, "The z component of the antiferromagnetic vector");
  params.addRequiredParam<Real>("beta", "The constant of spinflexoelectric interaction");
  params.addRequiredParam<Real>("P0", "The constant of the polarization scaling");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

SpinFlexoelectricPolarDerivative::SpinFlexoelectricPolarDerivative(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _polar_x_var(coupled("polar_x")),
  _polar_y_var(coupled("polar_y")),
  _polar_z_var(coupled("polar_z")),
  _antiferromag_L_x_var(coupled("antiferromag_L_x")),
  _antiferromag_L_y_var(coupled("antiferromag_L_y")),
  _antiferromag_L_z_var(coupled("antiferromag_L_z")),
  _antiferromag_L_x(coupledValue("antiferromag_L_x")),
  _antiferromag_L_y(coupledValue("antiferromag_L_y")),
  _antiferromag_L_z(coupledValue("antiferromag_L_z")),
  _antiferromag_L_x_grad(coupledGradient("antiferromag_L_x")),
  _antiferromag_L_y_grad(coupledGradient("antiferromag_L_y")),
  _antiferromag_L_z_grad(coupledGradient("antiferromag_L_z")),
  _beta(getParam<Real>("beta")),
  _P0(getParam<Real>("P0")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
SpinFlexoelectricPolarDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (_beta / _P0 ) * ((2*_antiferromag_L_x_grad[_qp](0) + _antiferromag_L_y_grad[_qp](1) + _antiferromag_L_z_grad[_qp](2))*_antiferromag_L_x[_qp] + _antiferromag_L_x_grad[_qp](1)*_antiferromag_L_y[_qp] + _antiferromag_L_x_grad[_qp](2)*_antiferromag_L_z[_qp]);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (_beta / _P0 ) * (_antiferromag_L_y_grad[_qp](0)*_antiferromag_L_x[_qp] + (_antiferromag_L_x_grad[_qp](0) + 2*_antiferromag_L_y_grad[_qp](1) + _antiferromag_L_z_grad[_qp](2))*_antiferromag_L_y[_qp] + _antiferromag_L_y_grad[_qp](2)*_antiferromag_L_z[_qp]);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (_beta / _P0 ) * (_antiferromag_L_z_grad[_qp](0)*_antiferromag_L_x[_qp] + _antiferromag_L_z_grad[_qp](1)*_antiferromag_L_y[_qp] + (_antiferromag_L_x_grad[_qp](0) + _antiferromag_L_y_grad[_qp](1) + 2*_antiferromag_L_z_grad[_qp](2))*_antiferromag_L_z[_qp]);
  }
  else
    return 0.0;
}

Real
SpinFlexoelectricPolarDerivative::computeQpJacobian()
{
  // no on-diagonal terms because the coupling is linear is P
  return 0.0;
}

Real
SpinFlexoelectricPolarDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiferromag_L_x_var)
    {
      return _test[_i][_qp] * (_beta / _P0 ) * ((2*_grad_phi[_j][_qp](0))*_antiferromag_L_x[_qp] + (2*_antiferromag_L_x_grad[_qp](0) + _antiferromag_L_y_grad[_qp](1) + _antiferromag_L_z_grad[_qp](2))*_phi[_j][_qp] + _grad_phi[_j][_qp](1)*_antiferromag_L_y[_qp] + _grad_phi[_j][_qp](2)*_antiferromag_L_z[_qp]);
    }
    else if (jvar == _antiferromag_L_y_var)
    {
      return _test[_i][_qp] * (_beta / _P0 ) * ((_grad_phi[_j][_qp](1))*_antiferromag_L_x[_qp] + _antiferromag_L_x_grad[_qp](1)*_phi[_j][_qp]);
    }
    else if (jvar == _antiferromag_L_z_var)
    {
      return _test[_i][_qp] * (_beta / _P0 ) * ((_grad_phi[_j][_qp](2))*_antiferromag_L_x[_qp] + _antiferromag_L_x_grad[_qp](2)*_phi[_j][_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _antiferromag_L_x_var)
    {
      return _test[_i][_qp] * (_beta / _P0 ) * (_antiferromag_L_y_grad[_qp](0)*_phi[_j][_qp] + (_grad_phi[_j][_qp](0))*_antiferromag_L_y[_qp]);
    }
    else if (jvar == _antiferromag_L_y_var)
    {
      return _test[_i][_qp] * (_beta / _P0 ) * (_grad_phi[_j][_qp](0)*_antiferromag_L_x[_qp] + (2*_grad_phi[_j][_qp](1))*_antiferromag_L_y[_qp] + (_antiferromag_L_x_grad[_qp](0) + 2*_antiferromag_L_y_grad[_qp](1) + _antiferromag_L_z_grad[_qp](2))*_phi[_j][_qp] + _grad_phi[_j][_qp](2)*_antiferromag_L_z[_qp]);
    }
    else if (jvar == _antiferromag_L_z_var)
    {
      return _test[_i][_qp] * (_beta / _P0 ) * ((_grad_phi[_j][_qp](2))*_antiferromag_L_y[_qp] + _antiferromag_L_y_grad[_qp](2)*_phi[_j][_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _antiferromag_L_x_var)
    {
      return _test[_i][_qp] * (_beta / _P0 ) * (_antiferromag_L_z_grad[_qp](0)*_phi[_j][_qp] + (_grad_phi[_j][_qp](0))*_antiferromag_L_z[_qp]);
    }
    else if (jvar == _antiferromag_L_y_var)
    {
      return _test[_i][_qp] * (_beta / _P0 ) * (_antiferromag_L_z_grad[_qp](1)*_phi[_j][_qp] + (_grad_phi[_j][_qp](1))*_antiferromag_L_z[_qp]);
    }
    else if (jvar == _antiferromag_L_z_var)
    {
      return _test[_i][_qp] * (_beta / _P0 ) * (_grad_phi[_j][_qp](0)*_antiferromag_L_x[_qp] + _grad_phi[_j][_qp](1)*_antiferromag_L_y[_qp] + (2*_grad_phi[_j][_qp](2))*_antiferromag_L_z[_qp] + (_antiferromag_L_x_grad[_qp](0) + _antiferromag_L_y_grad[_qp](1) + 2*_antiferromag_L_z_grad[_qp](2))*_phi[_j][_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
