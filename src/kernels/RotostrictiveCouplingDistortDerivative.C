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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "RotostrictiveCouplingDistortDerivative.h"

class RotostrictiveCouplingDistortDerivative;

registerMooseObject("FerretApp", RotostrictiveCouplingDistortDerivative);

template<>
InputParameters validParams<RotostrictiveCouplingDistortDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("u_x", "The x component of the local elastic displacement");
  params.addRequiredCoupledVar("u_y", "The y component of the local elastic displacement");
  params.addCoupledVar("u_z", 0.0, "The z component of the local elastic displacement");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the afd vector field");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the afd vector field");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the afd vector field");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

RotostrictiveCouplingDistortDerivative::RotostrictiveCouplingDistortDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _u_x_var(coupled("u_x")),
   _u_y_var(coupled("u_y")),
   _u_z_var(coupled("u_z")),
   _antiferrodis_A_x_var(coupled("antiferrodis_A_x")),
   _antiferrodis_A_y_var(coupled("antiferrodis_A_y")),
   _antiferrodis_A_z_var(coupled("antiferrodis_A_z")),
   _u_x_grad(coupledGradient("u_x")),
   _u_y_grad(coupledGradient("u_y")),
   _u_z_grad(coupledGradient("u_z")),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
   _r11(getMaterialProperty<Real>("r11")),
   _r12(getMaterialProperty<Real>("r12")),
   _r44(getMaterialProperty<Real>("r44"))
{
}

Real
RotostrictiveCouplingDistortDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return -0.5*_test[_i][_qp] * (-2*_antiferrodis_A_x[_qp]*_r11[_qp]*_u_x_grad[_qp](0) - 2*_r44[_qp]*((_antiferrodis_A_y[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0)))/2. + (_antiferrodis_A_z[_qp]*(_u_x_grad[_qp](2) + _u_z_grad[_qp](0)))/2.) - _r12[_qp]*(2*_antiferrodis_A_x[_qp]*_u_y_grad[_qp](1) + 2*_antiferrodis_A_x[_qp]*_u_z_grad[_qp](2)));
  }
  else if (_component == 1)
  {
    return -0.5*_test[_i][_qp] * (-2*_antiferrodis_A_y[_qp]*_r11[_qp]*_u_y_grad[_qp](1) - 2*_r44[_qp]*((_antiferrodis_A_x[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0)))/2. + (_antiferrodis_A_z[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)))/2.) - _r12[_qp]*(2*_antiferrodis_A_y[_qp]*_u_x_grad[_qp](0) + 2*_antiferrodis_A_y[_qp]*_u_z_grad[_qp](2)));
  }
  else if (_component == 2)
  {
    return -0.5*_test[_i][_qp] * (-(_r12[_qp]*(2*_antiferrodis_A_z[_qp]*_u_x_grad[_qp](0) + 2*_antiferrodis_A_z[_qp]*_u_y_grad[_qp](1))) - 2*_r44[_qp]*((_antiferrodis_A_x[_qp]*(_u_x_grad[_qp](2) + _u_z_grad[_qp](0)))/2. + (_antiferrodis_A_y[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)))/2.) - 2*_antiferrodis_A_z[_qp]*_r11[_qp]*_u_z_grad[_qp](2));
  }
  else
    return 0.0;
}

Real
RotostrictiveCouplingDistortDerivative::computeQpJacobian()
{
  if (_component == 0)
  {
    return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (-2*_r11[_qp]*_u_x_grad[_qp](0) - _r12[_qp]*(2*_u_y_grad[_qp](1) + 2*_u_z_grad[_qp](2)));
  }
  else if (_component == 1)
  {
    return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (-2*_r11[_qp]*_u_y_grad[_qp](1) - _r12[_qp]*(2*_u_x_grad[_qp](0) + 2*_u_z_grad[_qp](2)));
  }
  else if (_component == 2)
  {
    return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (-(_r12[_qp]*(2*_u_x_grad[_qp](0) + 2*_u_y_grad[_qp](1))) - 2*_r11[_qp]*_u_z_grad[_qp](2));
  }
  else
    return 0.0;
}

Real
RotostrictiveCouplingDistortDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiferrodis_A_y_var)
    {
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (-(_r44[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0))));
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (-(_r44[_qp]*(_u_x_grad[_qp](2) + _u_z_grad[_qp](0))));
    }
    else if (jvar == _u_x_var)
    {
      return -0.5*_test[_i][_qp] * (-2*_antiferrodis_A_x[_qp]*_grad_phi[_j][_qp](0)*_r11[_qp] - 2*((_antiferrodis_A_y[_qp]*_grad_phi[_j][_qp](1))/2. + (_antiferrodis_A_z[_qp]*_grad_phi[_j][_qp](2))/2.)*_r44[_qp]);
    }
    else if (jvar == _u_y_var)
    {
      return -0.5*_test[_i][_qp] * (-2*_antiferrodis_A_x[_qp]*_grad_phi[_j][_qp](1)*_r12[_qp] - _antiferrodis_A_y[_qp]*_grad_phi[_j][_qp](0)*_r44[_qp]);
    }
    else if (jvar == _u_z_var)
    {
      return -0.5*_test[_i][_qp] * (-2*_antiferrodis_A_x[_qp]*_grad_phi[_j][_qp](2)*_r12[_qp] - _antiferrodis_A_z[_qp]*_grad_phi[_j][_qp](0)*_r44[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return -0.5*_test[_i][_qp] * (-(_r44[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0))));
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return -0.5*_test[_i][_qp] * (-(_r44[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))));
    }
    else if (jvar == _u_x_var)
    {
      return -0.5*_test[_i][_qp] * (-2*_antiferrodis_A_y[_qp]*_grad_phi[_j][_qp](0)*_r12[_qp] - _antiferrodis_A_x[_qp]*_grad_phi[_j][_qp](1)*_r44[_qp]);
    }
    else if (jvar == _u_y_var)
    {
      return -0.5*_test[_i][_qp] * (-2*_antiferrodis_A_y[_qp]*_grad_phi[_j][_qp](1)*_r11[_qp] - 2*((_antiferrodis_A_x[_qp]*_grad_phi[_j][_qp](0))/2. + (_antiferrodis_A_z[_qp]*_grad_phi[_j][_qp](2))/2.)*_r44[_qp]);
    }
    else if (jvar == _u_z_var)
    {
      return -0.5*_test[_i][_qp] * (-2*_antiferrodis_A_y[_qp]*_grad_phi[_j][_qp](2)*_r12[_qp] - _antiferrodis_A_z[_qp]*_grad_phi[_j][_qp](1)*_r44[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return -0.5*_test[_i][_qp] * (-(_r44[_qp]*(_u_x_grad[_qp](2) + _u_z_grad[_qp](0))));
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return -0.5*_test[_i][_qp] * (-(_r44[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0))));
    }
    else if (jvar == _u_x_var)
    {
      return -0.5*_test[_i][_qp] * (-2*_antiferrodis_A_z[_qp]*_grad_phi[_j][_qp](0)*_r12[_qp] - _antiferrodis_A_x[_qp]*_grad_phi[_j][_qp](2)*_r44[_qp]);
    }
    else if (jvar == _u_y_var)
    {
      return -0.5*_test[_i][_qp] * (-2*_antiferrodis_A_z[_qp]*_grad_phi[_j][_qp](1)*_r12[_qp] - _antiferrodis_A_y[_qp]*_grad_phi[_j][_qp](2)*_r44[_qp]);
    }
    else if (jvar == _u_z_var)
    {
      return -0.5*_test[_i][_qp] * (-2*_antiferrodis_A_z[_qp]*_grad_phi[_j][_qp](2)*_r11[_qp] - 2*((_antiferrodis_A_x[_qp]*_grad_phi[_j][_qp](0))/2. + (_antiferrodis_A_y[_qp]*_grad_phi[_j][_qp](1))/2.)*_r44[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
