/**
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR _A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "DampedExchangeLLG.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<DampedExchangeLLG>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for the magnetic exchange coupling.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetic vector");
  params.addCoupledVar("mag_z", 0.0, "The z component of the magnetic vector");
  params.addRequiredParam<Real>("alphaLL", "the damping coefficient in the LLG equation");
  params.addRequiredParam<Real>("A", "exchange constant");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

DampedExchangeLLG::DampedExchangeLLG(const InputParameters & parameters)
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
  _alphaLL(getParam<Real>("alphaLL")),
  _A(getParam<Real>("A")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
DampedExchangeLLG::computeQpResidual()
{
  if (_component == 0)
  {
    return (2.0*_A*_alphaLL*_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp]*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp]*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp]*_mag_y[_qp] - 2.0*_A*_alphaLL*_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mag_y[_qp]) - 2.0*_A*_alphaLL*_mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mag_y[_qp]) - 
   2.0*_A*_alphaLL*_mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_y[_qp]) + 2.0*_A*_alphaLL*_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp]*_mag_z[_qp] - 2.0*_A*_alphaLL*_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mag_z[_qp]) - 
   2.0*_A*_alphaLL*_mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mag_z[_qp]) - 2.0*_A*_alphaLL*_mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_z[_qp]));
  }
  else if (_component == 1)
  {
    return (-2.0*_A*_alphaLL*_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mag_x[_qp]) - 2.0*_A*_alphaLL*_mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mag_x[_qp]) - 2.0*_A*_alphaLL*_mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_x[_qp]) + 2.0*_A*_alphaLL*_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp]*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp]*_mag_y[_qp] + 
   2.0*_A*_alphaLL*_mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp]*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp]*_mag_z[_qp] - 2.0*_A*_alphaLL*_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mag_z[_qp]) - 
   2.0*_A*_alphaLL*_mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mag_z[_qp]) - 2.0*_A*_alphaLL*_mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_z[_qp]));
  }
  else if (_component == 2)
  {
    return (-2.0*_A*_alphaLL*_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mag_x[_qp]) - 2.0*_A*_alphaLL*_mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mag_x[_qp]) - 2.0*_A*_alphaLL*_mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_x[_qp]) - 2.0*_A*_alphaLL*_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mag_y[_qp]) - 2.0*_A*_alphaLL*_mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mag_y[_qp]) - 2.0*_A*_alphaLL*_mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_y[_qp]) + 2.0*_A*_alphaLL*_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp]*_mag_z[_qp] + 
   2.0*_A*_alphaLL*_mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp]*_mag_z[_qp]);
  }
  else
    return 0.0;
}

Real
DampedExchangeLLG::computeQpJacobian()
{
  if (_component == 0)
  {
    return (2.0*_A*_alphaLL*_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_z[_qp] +(-2.0*_A*_alphaLL*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mag_y[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_y[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mag_z[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mag_z[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_z[_qp])));
  }
  else if (_component == 1)
  {
    return (2.0*_A*_alphaLL*_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_z[_qp]+(-2.0*_A*_alphaLL*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mag_x[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mag_x[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_x[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mag_z[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mag_z[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_z[_qp])));
  }
  else if (_component == 2)
  {
    return (2.0*_A*_alphaLL*_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp]+(-2.0*_A*_alphaLL*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mag_x[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mag_x[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_x[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*Utility::pow<2>(_mag_y[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*Utility::pow<2>(_mag_y[_qp]) - 2.0*_A*_alphaLL*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*Utility::pow<2>(_mag_y[_qp])));
  }
  else
    return 0.0;
}

Real
DampedExchangeLLG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return (2.0*_A*_alphaLL*_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp] - 4*_A*_alphaLL*_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp] - 4*_A*_alphaLL*_mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp] - 4*_A*_alphaLL*_mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp]+2.0*_A*_alphaLL*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp]*_mag_y[_qp] + 2.0*_A*_alphaLL*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp]*_mag_y[_qp] + 2.0*_A*_alphaLL*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp]*_mag_y[_qp]);
    }
    else if (jvar == _mag_z_var)
    {
      return (2.0*_A*_alphaLL*_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp] - 4*_A*_alphaLL*_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_z[_qp] - 4*_A*_alphaLL*_mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_z[_qp] - 4*_A*_alphaLL*_mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_z[_qp]+2.0*_A*_alphaLL*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp]*_mag_z[_qp]);
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
      return (-4*_A*_alphaLL*_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp] - 4*_A*_alphaLL*_mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp] - 4*_A*_alphaLL*_mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp]+2.0*_A*_alphaLL*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp]*_mag_y[_qp] + 2.0*_A*_alphaLL*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp]*_mag_y[_qp] + 2.0*_A*_alphaLL*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp]*_mag_y[_qp]);
    }
    else if (jvar == _mag_z_var)
    {
      return (2.0*_A*_alphaLL*_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp] - 4*_A*_alphaLL*_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_z[_qp] - 4*_A*_alphaLL*_mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_z[_qp] - 4*_A*_alphaLL*_mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_z[_qp]+2.0*_A*_alphaLL*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp]*_mag_z[_qp]);
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
      return (-4*_A*_alphaLL*_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp] - 4*_A*_alphaLL*_mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp] - 4*_A*_alphaLL*_mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_x_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_z[_qp]+2.0*_A*_alphaLL*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*_mag_x[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*_mag_x[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*_mag_x[_qp]*_mag_z[_qp]);
    }
    else if (jvar == _mag_y_var)
    {
      return (-4*_A*_alphaLL*_mag_z_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp] - 4*_A*_alphaLL*_mag_z_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp] - 4*_A*_alphaLL*_mag_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](0)*_grad_test[_i][_qp](0)*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_mag_z[_qp] + 2.0*_A*_alphaLL*_mag_y_grad[_qp](2)*_grad_test[_i][_qp](2)*_mag_z[_qp]+2.0*_A*_alphaLL*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*_mag_y[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*_mag_y[_qp]*_mag_z[_qp] + 2.0*_A*_alphaLL*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*_mag_y[_qp]*_mag_z[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
