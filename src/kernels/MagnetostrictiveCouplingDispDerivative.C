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

#include "MagnetostrictiveCouplingDispDerivative.h"
#include "libmesh/utility.h"

class MagnetostrictiveCouplingHeff;

registerMooseObject("FerretApp", MagnetostrictiveCouplingDispDerivative);

template<>
InputParameters validParams<MagnetostrictiveCouplingDispDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to the differentiation w.r.t spartial coordinates of the magnetoelastic self-strain"
                             " in the condition for mechanical equilibrium. Note for cubic magnets only.");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the mangetization");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredParam<Real>("C11", "the 11 (Voight) component of elastic stiffness tensor");
  params.addRequiredParam<Real>("C12", "the 12 (Voight) component of elastic stiffness tensor");
  params.addRequiredParam<Real>("C44", "the 44 (Voight) component of elastic stiffness tensor");
  params.addRequiredParam<Real>("L11", "the 11 (Voight) component of magnetostrictive coupling coefficient");
  params.addRequiredParam<Real>("L12", "the 12 (Voight) component of magnetostrictive coupling coefficient");
  params.addRequiredParam<Real>("L44", "the 44 (Voight) component of magnetostrictive coupling coefficient");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

MagnetostrictiveCouplingDispDerivative::MagnetostrictiveCouplingDispDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _mag_x_var(coupled("mag_x")),
   _mag_y_var(coupled("mag_y")),
   _mag_z_var(coupled("mag_z")),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _C11(getParam<Real>("C11")),
   _C12(getParam<Real>("C12")),
   _C44(getParam<Real>("C44")),
   _L11(getParam<Real>("L11")),
   _L12(getParam<Real>("L12")),
   _L44(getParam<Real>("L44"))
{
}

Real
MagnetostrictiveCouplingDispDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return ((_C11*Utility::pow<2>(_mag_x[_qp])*_L11 + _C12*Utility::pow<2>(_mag_y[_qp])*_L11 + _C12*Utility::pow<2>(_mag_z[_qp])*_L11 + 
      2*_C12*Utility::pow<2>(_mag_x[_qp])*_L12 + _C11*Utility::pow<2>(_mag_y[_qp])*_L12 + _C12*Utility::pow<2>(_mag_y[_qp])*_L12 + 
      _C11*Utility::pow<2>(_mag_z[_qp])*_L12 + _C12*Utility::pow<2>(_mag_z[_qp])*_L12)*_grad_test[_i][_qp](0) + 4*_C44*_mag_x[_qp]*_mag_y[_qp]*_L44*_grad_test[_i][_qp](1) + 
   4*_C44*_mag_x[_qp]*_mag_z[_qp]*_L44*_grad_test[_i][_qp](2));
  }
  else if (_component == 1)
  {
    return (4*_C44*_mag_x[_qp]*_mag_y[_qp]*_L44*_grad_test[_i][_qp](0) + (_C12*Utility::pow<2>(_mag_x[_qp])*_L11 + _C11*Utility::pow<2>(_mag_y[_qp])*_L11 + 
      _C12*Utility::pow<2>(_mag_z[_qp])*_L11 + _C11*Utility::pow<2>(_mag_x[_qp])*_L12 + _C12*Utility::pow<2>(_mag_x[_qp])*_L12 + 
      2*_C12*Utility::pow<2>(_mag_y[_qp])*_L12 + _C11*Utility::pow<2>(_mag_z[_qp])*_L12 + _C12*Utility::pow<2>(_mag_z[_qp])*_L12)*_grad_test[_i][_qp](1) + 
   4*_C44*_mag_y[_qp]*_mag_z[_qp]*_L44*_grad_test[_i][_qp](2));
  }
  else if (_component == 2)
  {
    return (4*_C44*_mag_x[_qp]*_mag_z[_qp]*_L44*_grad_test[_i][_qp](0) + 4*_C44*_mag_y[_qp]*_mag_z[_qp]*_L44*_grad_test[_i][_qp](1) + 
   (_C12*Utility::pow<2>(_mag_x[_qp])*_L11 + _C12*Utility::pow<2>(_mag_y[_qp])*_L11 + _C11*Utility::pow<2>(_mag_z[_qp])*_L11 + 
      _C11*Utility::pow<2>(_mag_x[_qp])*_L12 + _C12*Utility::pow<2>(_mag_x[_qp])*_L12 + _C11*Utility::pow<2>(_mag_y[_qp])*_L12 + 
      _C12*Utility::pow<2>(_mag_y[_qp])*_L12 + 2*_C12*Utility::pow<2>(_mag_z[_qp])*_L12)*_grad_test[_i][_qp](2));
  }
  else
    return 0.0;
}

Real
MagnetostrictiveCouplingDispDerivative::computeQpJacobian()
{
  return 0.0;
}

Real
MagnetostrictiveCouplingDispDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_x_var)
    {
      return _phi[_j][_qp] * (2*_C11*_mag_x[_qp]*_L11*_grad_test[_i][_qp](0) + 4*_C12*_mag_x[_qp]*_L12*_grad_test[_i][_qp](0) + 4*_C44*_mag_y[_qp]*_L44*_grad_test[_i][_qp](1) + 4*_C44*_mag_z[_qp]*_L44*_grad_test[_i][_qp](2));
    }
    else if (jvar == _mag_y_var)
    {
      return _phi[_j][_qp] * (2*_C12*_mag_y[_qp]*_L11*_grad_test[_i][_qp](0) + 2*_C11*_mag_y[_qp]*_L12*_grad_test[_i][_qp](0) + 2*_C12*_mag_y[_qp]*_L12*_grad_test[_i][_qp](0) + 4*_C44*_mag_x[_qp]*_L44*_grad_test[_i][_qp](1));
    }
    else if (jvar == _mag_z_var)
    {
      return _phi[_j][_qp] * (2*_C12*_mag_z[_qp]*_L11*_grad_test[_i][_qp](0) + 2*_C11*_mag_z[_qp]*_L12*_grad_test[_i][_qp](0) + 2*_C12*_mag_z[_qp]*_L12*_grad_test[_i][_qp](0) + 4*_C44*_mag_x[_qp]*_L44*_grad_test[_i][_qp](2));
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
      return _phi[_j][_qp] * (4*_C44*_mag_y[_qp]*_L44*_grad_test[_i][_qp](0) + 2*_C12*_mag_x[_qp]*_L11*_grad_test[_i][_qp](1) + 2*_C11*_mag_x[_qp]*_L12*_grad_test[_i][_qp](1) + 2*_C12*_mag_x[_qp]*_L12*_grad_test[_i][_qp](1));
    }
    else if (jvar == _mag_y_var)
    {
      return _phi[_j][_qp] * (4*_C44*_mag_x[_qp]*_L44*_grad_test[_i][_qp](0) + 2*_C11*_mag_y[_qp]*_L11*_grad_test[_i][_qp](1) + 4*_C12*_mag_y[_qp]*_L12*_grad_test[_i][_qp](1) + 4*_C44*_mag_z[_qp]*_L44*_grad_test[_i][_qp](2));
    }
    else if (jvar == _mag_z_var)
    {
      return _phi[_j][_qp] * (2*_C12*_mag_z[_qp]*_L11*_grad_test[_i][_qp](1) + 2*_C11*_mag_z[_qp]*_L12*_grad_test[_i][_qp](1) + 2*_C12*_mag_z[_qp]*_L12*_grad_test[_i][_qp](1) + 4*_C44*_mag_y[_qp]*_L44*_grad_test[_i][_qp](2));
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
      return _phi[_j][_qp] * (4*_C44*_mag_z[_qp]*_L44*_grad_test[_i][_qp](0) + 2*_C12*_mag_x[_qp]*_L11*_grad_test[_i][_qp](2) + 2*_C11*_mag_x[_qp]*_L12*_grad_test[_i][_qp](2) + 2*_C12*_mag_x[_qp]*_L12*_grad_test[_i][_qp](2));
    }
    else if (jvar == _mag_y_var)
    {
      return _phi[_j][_qp] * (4*_C44*_mag_z[_qp]*_L44*_grad_test[_i][_qp](1) + 2*_C12*_mag_y[_qp]*_L11*_grad_test[_i][_qp](2) + 2*_C11*_mag_y[_qp]*_L12*_grad_test[_i][_qp](2) + 2*_C12*_mag_y[_qp]*_L12*_grad_test[_i][_qp](2));
    }
    else if (jvar == _mag_z_var)
    {
      return _phi[_j][_qp] * (4*_C44*_mag_x[_qp]*_L44*_grad_test[_i][_qp](0) + 4*_C44*_mag_y[_qp]*_L44*_grad_test[_i][_qp](1) + 2*_C11*_mag_z[_qp]*_L11*_grad_test[_i][_qp](2) + 4*_C12*_mag_z[_qp]*_L12*_grad_test[_i][_qp](2));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
