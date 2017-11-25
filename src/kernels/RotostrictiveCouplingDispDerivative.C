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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "RotostrictiveCouplingDispDerivative.h"

class RotostrictiveCouplingDispDerivative;

template<>
InputParameters validParams<RotostrictiveCouplingDispDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the afd vector field");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the afd vector field");
  params.addRequiredCoupledVar("antiferrodis_A_z", "The z component of the afd vector field");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredParam<Real>("C11", "the 11 (Voight) component of elastic stiffness tensor");
  params.addRequiredParam<Real>("C12", "the 12 (Voight) component of elastic stiffness tensor");
  params.addRequiredParam<Real>("C44", "the 44 (Voight) component of elastic stiffness tensor");
  params.addRequiredParam<Real>("R11", "the 11 (Voight) component of rotostrictive coupling coefficient");
  params.addRequiredParam<Real>("R12", "the 12 (Voight) component of rotostrictive coupling coefficient");
  params.addRequiredParam<Real>("R44", "the 44 (Voight) component of rotostrictive coupling coefficient");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

RotostrictiveCouplingDispDerivative::RotostrictiveCouplingDispDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _antiferrodis_A_x_var(coupled("antiferrodis_A_x")),
   _antiferrodis_A_y_var(coupled("antiferrodis_A_y")),
   _antiferrodis_A_z_var(coupled("antiferrodis_A_z")),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
   _C11(getParam<Real>("C11")),
   _C12(getParam<Real>("C12")),
   _C44(getParam<Real>("C44")),
   _R11(getParam<Real>("R11")),
   _R12(getParam<Real>("R12")),
   _R44(getParam<Real>("R44")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
RotostrictiveCouplingDispDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return -(_grad_test[_i][_qp](0)*(_C12*(std::pow(_antiferrodis_A_z[_qp],2)*_R11 + (std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*_R12) + _C12*(std::pow(_antiferrodis_A_y[_qp],2)*_R11 + (std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_R12) + _C11*(std::pow(_antiferrodis_A_x[_qp],2)*_R11 + (std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_R12)) + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_C44*_grad_test[_i][_qp](1)*_R44 + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_C44*_grad_test[_i][_qp](2)*_R44);
  }
  else if (_component == 1)
  {
    return -(_grad_test[_i][_qp](1)*(_C12*(std::pow(_antiferrodis_A_z[_qp],2)*_R11 + (std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*_R12) + _C11*(std::pow(_antiferrodis_A_y[_qp],2)*_R11 + (std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_R12) + _C12*(std::pow(_antiferrodis_A_x[_qp],2)*_R11 + (std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_R12)) + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_C44*_grad_test[_i][_qp](0)*_R44 + 2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_C44*_grad_test[_i][_qp](2)*_R44);
  }
  else if (_component == 2)
  {
    return -(_grad_test[_i][_qp](2)*(_C11*(std::pow(_antiferrodis_A_z[_qp],2)*_R11 + (std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*_R12) + _C12*(std::pow(_antiferrodis_A_y[_qp],2)*_R11 + (std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_R12) + _C12*(std::pow(_antiferrodis_A_x[_qp],2)*_R11 + (std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_R12)) + 
   2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_C44*_grad_test[_i][_qp](0)*_R44 + 2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_C44*_grad_test[_i][_qp](1)*_R44);
  }
  else
    return 0.0;
}

Real
RotostrictiveCouplingDispDerivative::computeQpJacobian()
{
  return 0.0;
}

Real
RotostrictiveCouplingDispDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return -_phi[_j][_qp] * (_grad_test[_i][_qp](0)*(2*_antiferrodis_A_x[_qp]*_C11*_R11 + 4*_antiferrodis_A_x[_qp]*_C12*_R12) + 2*_antiferrodis_A_y[_qp]*_C44*_grad_test[_i][_qp](1)*_R44 + 2*_antiferrodis_A_z[_qp]*_C44*_grad_test[_i][_qp](2)*_R44);
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return -_phi[_j][_qp] * (_grad_test[_i][_qp](0)*(2*_antiferrodis_A_y[_qp]*_C12*_R11 + 2*_antiferrodis_A_y[_qp]*_C11*_R12 + 2*_antiferrodis_A_y[_qp]*_C12*_R12) + 2*_antiferrodis_A_x[_qp]*_C44*_grad_test[_i][_qp](1)*_R44);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return -_phi[_j][_qp] * (_grad_test[_i][_qp](0)*(2*_antiferrodis_A_z[_qp]*_C12*_R11 + 2*_antiferrodis_A_z[_qp]*_C11*_R12 + 2*_antiferrodis_A_z[_qp]*_C12*_R12) + 2*_antiferrodis_A_x[_qp]*_C44*_grad_test[_i][_qp](2)*_R44);
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
      return -_phi[_j][_qp] * (_grad_test[_i][_qp](1)*(2*_antiferrodis_A_x[_qp]*_C12*_R11 + 2*_antiferrodis_A_x[_qp]*_C11*_R12 + 2*_antiferrodis_A_x[_qp]*_C12*_R12) + 2*_antiferrodis_A_y[_qp]*_C44*_grad_test[_i][_qp](0)*_R44);
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return -_phi[_j][_qp] * (_grad_test[_i][_qp](1)*(2*_antiferrodis_A_y[_qp]*_C11*_R11 + 4*_antiferrodis_A_y[_qp]*_C12*_R12) + 2*_antiferrodis_A_x[_qp]*_C44*_grad_test[_i][_qp](0)*_R44 + 2*_antiferrodis_A_z[_qp]*_C44*_grad_test[_i][_qp](2)*_R44);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return -_phi[_j][_qp] * (_grad_test[_i][_qp](1)*(2*_antiferrodis_A_z[_qp]*_C12*_R11 + 2*_antiferrodis_A_z[_qp]*_C11*_R12 + 2*_antiferrodis_A_z[_qp]*_C12*_R12) + 2*_antiferrodis_A_y[_qp]*_C44*_grad_test[_i][_qp](2)*_R44);
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
      return -_phi[_j][_qp] * (_grad_test[_i][_qp](2)*(2*_antiferrodis_A_x[_qp]*_C12*_R11 + 2*_antiferrodis_A_x[_qp]*_C11*_R12 + 2*_antiferrodis_A_x[_qp]*_C12*_R12) + 2*_antiferrodis_A_z[_qp]*_C44*_grad_test[_i][_qp](0)*_R44);
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return -_phi[_j][_qp] * (_grad_test[_i][_qp](2)*(2*_antiferrodis_A_y[_qp]*_C12*_R11 + 2*_antiferrodis_A_y[_qp]*_C11*_R12 + 2*_antiferrodis_A_y[_qp]*_C12*_R12) + 2*_antiferrodis_A_z[_qp]*_C44*_grad_test[_i][_qp](1)*_R44);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return -_phi[_j][_qp] * (_grad_test[_i][_qp](2)*(2*_antiferrodis_A_z[_qp]*_C11*_R11 + 4*_antiferrodis_A_z[_qp]*_C12*_R12) + 2*_antiferrodis_A_x[_qp]*_C44*_grad_test[_i][_qp](0)*_R44 + 2*_antiferrodis_A_y[_qp]*_C44*_grad_test[_i][_qp](1)*_R44);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
