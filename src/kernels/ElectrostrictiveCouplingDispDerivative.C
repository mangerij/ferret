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

#include "ElectrostrictiveCouplingDispDerivative.h"

class ElectrostrictiveCouplingDispDerivative;

template<>
InputParameters validParams<ElectrostrictiveCouplingDispDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredParam<Real>("C11", "the 11 (Voight) component of elastic stiffness tensor");
  params.addRequiredParam<Real>("C12", "the 12 (Voight) component of elastic stiffness tensor");
  params.addRequiredParam<Real>("C44", "the 44 (Voight) component of elastic stiffness tensor");
  params.addRequiredParam<Real>("Q11", "the 11 (Voight) component of electrostrictive coupling coefficient");
  params.addRequiredParam<Real>("Q12", "the 12 (Voight) component of electrostrictive coupling coefficient");
  params.addRequiredParam<Real>("Q44", "the 44 (Voight) component of electrostrictive coupling coefficient");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

ElectrostrictiveCouplingDispDerivative::ElectrostrictiveCouplingDispDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _C11(getParam<Real>("C11")),
   _C12(getParam<Real>("C12")),
   _C44(getParam<Real>("C44")),
   _Q11(getParam<Real>("Q11")),
   _Q12(getParam<Real>("Q12")),
   _Q44(getParam<Real>("Q44")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
ElectrostrictiveCouplingDispDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return (2*_polar_x[_qp]*_polar_y[_qp]*_C44*_grad_test[_i][_qp](1)*_Q44 + 2*_polar_x[_qp]*_polar_z[_qp]*_C44*_grad_test[_i][_qp](2)*_Q44 + _grad_test[_i][_qp](0)*(_C12*(std::pow(_polar_z[_qp],2)*_Q11 + (std::pow(_polar_x[_qp],2) + std::pow(_polar_y[_qp],2))*_Q12) + _C12*(std::pow(_polar_y[_qp],2)*_Q11 + (std::pow(_polar_x[_qp],2) + std::pow(_polar_z[_qp],2))*_Q12) + 
      _C11*((std::pow(_polar_y[_qp],2) + std::pow(_polar_z[_qp],2))*_Q12 + std::pow(_polar_x[_qp],2)*_Q44)));
  }
  else if (_component == 1)
  {
    return (2*_polar_x[_qp]*_polar_y[_qp]*_C44*_grad_test[_i][_qp](0)*_Q44 + 2*_polar_y[_qp]*_polar_z[_qp]*_C44*_grad_test[_i][_qp](2)*_Q44 + _grad_test[_i][_qp](1)*(_C12*(std::pow(_polar_z[_qp],2)*_Q11 + (std::pow(_polar_x[_qp],2) + std::pow(_polar_y[_qp],2))*_Q12) + _C11*(std::pow(_polar_y[_qp],2)*_Q11 + (std::pow(_polar_x[_qp],2) + std::pow(_polar_z[_qp],2))*_Q12) + 
      _C12*((std::pow(_polar_y[_qp],2) + std::pow(_polar_z[_qp],2))*_Q12 + std::pow(_polar_x[_qp],2)*_Q44)));
  }
  else if (_component == 2)
  {
    return (2*_polar_x[_qp]*_polar_z[_qp]*_C44*_grad_test[_i][_qp](0)*_Q44 + 2*_polar_y[_qp]*_polar_z[_qp]*_C44*_grad_test[_i][_qp](1)*_Q44 + _grad_test[_i][_qp](2)*(_C11*(std::pow(_polar_z[_qp],2)*_Q11 + (std::pow(_polar_x[_qp],2) + std::pow(_polar_y[_qp],2))*_Q12) + _C12*(std::pow(_polar_y[_qp],2)*_Q11 + (std::pow(_polar_x[_qp],2) + std::pow(_polar_z[_qp],2))*_Q12) + 
      _C12*((std::pow(_polar_y[_qp],2) + std::pow(_polar_z[_qp],2))*_Q12 + std::pow(_polar_x[_qp],2)*_Q44)));
  }
  else
    return 0.0;
}

Real
ElectrostrictiveCouplingDispDerivative::computeQpJacobian()
{
  return 0.0;
}

Real
ElectrostrictiveCouplingDispDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_x_var)
    {
      return _phi[_j][_qp] * (2*_polar_y[_qp]*_C44*_grad_test[_i][_qp](1)*_Q44 + 2*_polar_z[_qp]*_C44*_grad_test[_i][_qp](2)*_Q44 + _grad_test[_i][_qp](0)*(4*_polar_x[_qp]*_C12*_Q12 + 2*_polar_x[_qp]*_C11*_Q44));
    }
    else if (jvar == _polar_y_var)
    {
      return _phi[_j][_qp] * (_grad_test[_i][_qp](0)*(2*_polar_y[_qp]*_C12*_Q11 + 2*_polar_y[_qp]*_C11*_Q12 + 2*_polar_y[_qp]*_C12*_Q12) + 2*_polar_x[_qp]*_C44*_grad_test[_i][_qp](1)*_Q44);
    }
    else if (jvar == _polar_z_var)
    {
      return _phi[_j][_qp] * (_grad_test[_i][_qp](0)*(2*_polar_z[_qp]*_C12*_Q11 + 2*_polar_z[_qp]*_C11*_Q12 + 2*_polar_z[_qp]*_C12*_Q12) + 2*_polar_x[_qp]*_C44*_grad_test[_i][_qp](2)*_Q44);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _polar_x_var)
    {
      return _phi[_j][_qp] * (2*_polar_y[_qp]*_C44*_grad_test[_i][_qp](0)*_Q44 + _grad_test[_i][_qp](1)*(2*_polar_x[_qp]*_C11*_Q12 + 2*_polar_x[_qp]*_C12*_Q12 + 2*_polar_x[_qp]*_C12*_Q44));
    }
    else if (jvar == _polar_y_var)
    {
      return _phi[_j][_qp] * (_grad_test[_i][_qp](1)*(2*_polar_y[_qp]*_C11*_Q11 + 4*_polar_y[_qp]*_C12*_Q12) + 2*_polar_x[_qp]*_C44*_grad_test[_i][_qp](0)*_Q44 + 2*_polar_z[_qp]*_C44*_grad_test[_i][_qp](2)*_Q44);
    }
    else if (jvar == _polar_z_var)
    {
      return _phi[_j][_qp] * (_grad_test[_i][_qp](1)*(2*_polar_z[_qp]*_C12*_Q11 + 2*_polar_z[_qp]*_C11*_Q12 + 2*_polar_z[_qp]*_C12*_Q12) + 2*_polar_y[_qp]*_C44*_grad_test[_i][_qp](2)*_Q44);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _polar_x_var)
    {
      return _phi[_j][_qp] * (2*_polar_z[_qp]*_C44*_grad_test[_i][_qp](0)*_Q44 + _grad_test[_i][_qp](2)*(2*_polar_x[_qp]*_C11*_Q12 + 2*_polar_x[_qp]*_C12*_Q12 + 2*_polar_x[_qp]*_C12*_Q44));
    }
    else if (jvar == _polar_y_var)
    {
      return _phi[_j][_qp] * (_grad_test[_i][_qp](2)*(2*_polar_y[_qp]*_C12*_Q11 + 2*_polar_y[_qp]*_C11*_Q12 + 2*_polar_y[_qp]*_C12*_Q12) + 2*_polar_z[_qp]*_C44*_grad_test[_i][_qp](1)*_Q44);
    }
    else if (jvar == _polar_z_var)
    {
      return _phi[_j][_qp] * (_grad_test[_i][_qp](2)*(2*_polar_z[_qp]*_C11*_Q11 + 4*_polar_z[_qp]*_C12*_Q12) + 2*_polar_x[_qp]*_C44*_grad_test[_i][_qp](0)*_Q44 + 2*_polar_y[_qp]*_C44*_grad_test[_i][_qp](1)*_Q44);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
