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

#include "MagnetostrictiveCouplingCubicHeff.h"
#include "libmesh/utility.h"

class MagnetostrictiveCouplingCubicHeff;

registerMooseObject("FerretApp", MagnetostrictiveCouplingCubicHeff);

template<>
InputParameters validParams<MagnetostrictiveCouplingCubicHeff>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to the magnetoelectric effective field. Note for cubic magnets only.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the magnetization");
  params.addRequiredCoupledVar("disp_x", "The x component of the elastic displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the elastic displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the elastic displacement");
  return params;
}

MagnetostrictiveCouplingCubicHeff::MagnetostrictiveCouplingCubicHeff(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _mag_x_var(coupled("mag_x")),
   _mag_y_var(coupled("mag_y")),
   _mag_z_var(coupled("mag_z")),
   _disp_x_var(coupled("disp_x")),
   _disp_y_var(coupled("disp_y")),
   _disp_z_var(coupled("disp_z")),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _disp_x_grad(coupledGradient("disp_x")),
   _disp_y_grad(coupledGradient("disp_y")),
   _disp_z_grad(coupledGradient("disp_z")),
   _alpha(getMaterialProperty<Real>("_alpha")),
   _g0(getMaterialProperty<Real>("g0")),
   _Ms(getMaterialProperty<Real>("Ms")),
   _B11(getMaterialProperty<Real>("B11")),
   _B22(getMaterialProperty<Real>("B22"))
{
}

Real
MagnetostrictiveCouplingCubicHeff::computeQpResidual()
{
  if (_component == 0)
  {
    return (_g0[_qp]*(2*_B11[_qp]*((-_disp_y_grad[_qp](1) + _disp_z_grad[_qp](2))*_mag_y[_qp]*_mag_z[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*((_disp_x_grad[_qp](0) - _disp_y_grad[_qp](1))*Utility::pow<2>(_mag_y[_qp]) + (_disp_x_grad[_qp](0) - _disp_z_grad[_qp](2))*Utility::pow<2>(_mag_z[_qp]))) + 
       _B22[_qp]*(-(0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*_mag_x[_qp]*_mag_z[_qp]) + 0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*(_mag_y[_qp] - _mag_z[_qp])*(_mag_y[_qp] + _mag_z[_qp]) + 0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*(_mag_x[_qp]*_mag_y[_qp] - _alpha[_qp]*_Ms[_qp]*Utility::pow<2>(_mag_x[_qp])*_mag_z[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_z[_qp]*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))) + _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*(-2*0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_x[_qp]*_mag_z[_qp] + 0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*(-Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))))*
     _test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
  }
  else if (_component == 1)
  {
    return (_g0[_qp]*(_B22[_qp]*(-(0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_x[_qp]*_mag_y[_qp]) + 0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*_mag_y[_qp]*_mag_z[_qp] + _alpha[_qp]*_Ms[_qp]*(0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*_mag_x[_qp] + 0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_z[_qp])*(Utility::pow<2>(_mag_x[_qp]) - Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])) + 0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*(-Utility::pow<2>(_mag_x[_qp]) - 2*_alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_y[_qp]*_mag_z[_qp] + Utility::pow<2>(_mag_z[_qp]))) + 
       2*_B11[_qp]*((_disp_x_grad[_qp](0) - _disp_z_grad[_qp](2))*_mag_x[_qp]*_mag_z[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*((-_disp_x_grad[_qp](0) + _disp_y_grad[_qp](1))*Utility::pow<2>(_mag_x[_qp]) + (_disp_y_grad[_qp](1) - _disp_z_grad[_qp](2))*Utility::pow<2>(_mag_z[_qp]))))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
  }
  else if (_component == 2)
  {
    return (_g0[_qp]*(2*_B11[_qp]*((-_disp_x_grad[_qp](0) + _disp_y_grad[_qp](1))*_mag_x[_qp]*_mag_y[_qp] + _alpha[_qp]*_Ms[_qp]*((-_disp_x_grad[_qp](0) + _disp_z_grad[_qp](2))*Utility::pow<2>(_mag_x[_qp]) + (-_disp_y_grad[_qp](1) + _disp_z_grad[_qp](2))*Utility::pow<2>(_mag_y[_qp]))*_mag_z[_qp]) + 
       _B22[_qp]*((0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_x[_qp] - 0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*_mag_y[_qp])*_mag_z[_qp] + 0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*(Utility::pow<2>(_mag_x[_qp]) - Utility::pow<2>(_mag_y[_qp]) - 2*_alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_y[_qp]*_mag_z[_qp]) + _alpha[_qp]*_Ms[_qp]*(0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*_mag_x[_qp] + 0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_y[_qp])*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) - Utility::pow<2>(_mag_z[_qp]))))*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
  }
  else
    return 0.0;
}

Real
MagnetostrictiveCouplingCubicHeff::computeQpJacobian()
{
  if (_component == 0)
  {
    return (_g0[_qp]*(_B22[_qp]*(0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0)) - 2*_alpha[_qp]*0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*_Ms[_qp]*_mag_x[_qp])*_mag_y[_qp] - _B22[_qp]*(0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0)) + 2*_alpha[_qp]*_Ms[_qp]*(0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*_mag_x[_qp] + 0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_y[_qp]))*_mag_z[_qp] + 2*_alpha[_qp]*_B11[_qp]*_Ms[_qp]*((_disp_x_grad[_qp](0) - _disp_y_grad[_qp](1))*Utility::pow<2>(_mag_y[_qp]) + (_disp_x_grad[_qp](0) - _disp_z_grad[_qp](2))*Utility::pow<2>(_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
  }
  else if (_component == 1)
  {
    return -((_g0[_qp]*(2*_alpha[_qp]*_B11[_qp]*_Ms[_qp]*((_disp_x_grad[_qp](0) - _disp_y_grad[_qp](1))*Utility::pow<2>(_mag_x[_qp]) + (-_disp_y_grad[_qp](1) + _disp_z_grad[_qp](2))*Utility::pow<2>(_mag_z[_qp])) + _B22[_qp]*(-(0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*_mag_z[_qp]) + 2*_alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*(0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*_mag_y[_qp] + 0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*_mag_z[_qp]) + 0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*(_mag_x[_qp] + 2*_alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp])))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp])));
  }
  else if (_component == 2)
  {
    return -((_g0[_qp]*(2*_alpha[_qp]*_B11[_qp]*_Ms[_qp]*((_disp_x_grad[_qp](0) - _disp_z_grad[_qp](2))*Utility::pow<2>(_mag_x[_qp]) + (_disp_y_grad[_qp](1) - _disp_z_grad[_qp](2))*Utility::pow<2>(_mag_y[_qp])) + _B22[_qp]*(-(0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_x[_qp]) + 0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*_mag_y[_qp] + 2*_alpha[_qp]*0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*_Ms[_qp]*_mag_x[_qp]*_mag_y[_qp] + 2*_alpha[_qp]*_Ms[_qp]*(0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*_mag_x[_qp] + 0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_y[_qp])*_mag_z[_qp]))*_phi[_j][_qp]*_test[_i][_qp])/((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp])));
  }
  else
    return 0.0;
}

Real
MagnetostrictiveCouplingCubicHeff::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _mag_y_var)
    {
      return (_g0[_qp]*(4*_alpha[_qp]*_B11[_qp]*(_disp_x_grad[_qp](0) - _disp_y_grad[_qp](1))*_Ms[_qp]*_mag_x[_qp]*_mag_y[_qp] + 2*_B11[_qp]*(-_disp_y_grad[_qp](1) + _disp_z_grad[_qp](2))*_mag_z[_qp] + _B22[_qp]*(2*0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_y[_qp] + 0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*(_mag_x[_qp] + 2*_alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp]) + _alpha[_qp]*_Ms[_qp]*(-2*0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_x[_qp]*_mag_z[_qp] + 0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*(-Utility::pow<2>(_mag_x[_qp]) + 3*Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))))*_phi[_j][_qp]*_test[_i][_qp])/
   ((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
    }
    else if (jvar == _mag_z_var)
    {
      return -((_g0[_qp]*(2*_B11[_qp]*(_disp_y_grad[_qp](1) - _disp_z_grad[_qp](2))*_mag_y[_qp] + 4*_alpha[_qp]*_B11[_qp]*(-_disp_x_grad[_qp](0) + _disp_z_grad[_qp](2))*_Ms[_qp]*_mag_x[_qp]*_mag_z[_qp] + _B22[_qp]*(2*0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_z[_qp] + 0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*(_mag_x[_qp] - 2*_alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp]) + _alpha[_qp]*_Ms[_qp]*(2*0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_x[_qp]*_mag_y[_qp] + 0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*(Utility::pow<2>(_mag_x[_qp]) - Utility::pow<2>(_mag_y[_qp]) - 3*Utility::pow<2>(_mag_z[_qp])))))*_phi[_j][_qp]*_test[_i][_qp])/
     ((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp])));
    }
    else if (jvar == _disp_x_var)
    {
      return (_g0[_qp]*(4*_alpha[_qp]*_B11[_qp]*_grad_phi[_j][_qp](0)*_Ms[_qp]*_mag_x[_qp]*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])) + _B22[_qp]*(_grad_phi[_j][_qp](2)*(_mag_x[_qp]*_mag_y[_qp] - _alpha[_qp]*_Ms[_qp]*Utility::pow<2>(_mag_x[_qp])*_mag_z[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_z[_qp]*(Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))) + 
          _grad_phi[_j][_qp](1)*(-(_mag_x[_qp]*_mag_z[_qp]) + _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*(-Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))))*_test[_i][_qp])/(2.*(1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
    }
    else if (jvar == _disp_y_var)
    {
      return -(_g0[_qp]*(4*_B11[_qp]*_grad_phi[_j][_qp](1)*_mag_y[_qp]*_mag_z[_qp] + _B22[_qp]*(-(_grad_phi[_j][_qp](2)*Utility::pow<2>(_mag_y[_qp])) + _grad_phi[_j][_qp](0)*_mag_x[_qp]*_mag_z[_qp] + _grad_phi[_j][_qp](2)*Utility::pow<2>(_mag_z[_qp])) + _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*(4*_B11[_qp]*_grad_phi[_j][_qp](1)*_mag_x[_qp]*_mag_y[_qp] + _B22[_qp]*(2*_grad_phi[_j][_qp](2)*_mag_x[_qp]*_mag_z[_qp] + _grad_phi[_j][_qp](0)*(Utility::pow<2>(_mag_x[_qp]) - Utility::pow<2>(_mag_y[_qp]) - Utility::pow<2>(_mag_z[_qp])))))*_test[_i][_qp])/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
    }
    else if (jvar == _disp_z_var)
    {
      return (_g0[_qp]*(4*_B11[_qp]*_grad_phi[_j][_qp](2)*_mag_z[_qp]*(_mag_y[_qp] - _alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_z[_qp]) + _B22[_qp]*(_mag_y[_qp]*(_grad_phi[_j][_qp](0)*_mag_x[_qp] + _grad_phi[_j][_qp](1)*_mag_y[_qp]) + _alpha[_qp]*_Ms[_qp]*(-2*_grad_phi[_j][_qp](1)*_mag_x[_qp]*_mag_y[_qp] + _grad_phi[_j][_qp](0)*(-Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp])))*_mag_z[_qp] - _grad_phi[_j][_qp](1)*Utility::pow<2>(_mag_z[_qp]) + _alpha[_qp]*_grad_phi[_j][_qp](0)*_Ms[_qp]*Utility::pow<3>(_mag_z[_qp])))*_test[_i][_qp])/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
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
      return -((_g0[_qp]*(4*_alpha[_qp]*_B11[_qp]*(_disp_x_grad[_qp](0) - _disp_y_grad[_qp](1))*_Ms[_qp]*_mag_x[_qp]*_mag_y[_qp] + 2*_B11[_qp]*(-_disp_x_grad[_qp](0) + _disp_z_grad[_qp](2))*_mag_z[_qp] + _B22[_qp]*(0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_y[_qp] + 2*0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*(_mag_x[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp]) - _alpha[_qp]*_Ms[_qp]*(2*0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_x[_qp]*_mag_z[_qp] + 0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*(3*Utility::pow<2>(_mag_x[_qp]) - Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))))*_phi[_j][_qp]*_test[_i][_qp])/
     ((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp])));
    }
    else if (jvar == _mag_z_var)
    {
      return (_g0[_qp]*(2*_B11[_qp]*(_disp_x_grad[_qp](0) - _disp_z_grad[_qp](2))*_mag_x[_qp] + 4*_alpha[_qp]*_B11[_qp]*(_disp_y_grad[_qp](1) - _disp_z_grad[_qp](2))*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp] + _B22[_qp]*(0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*_mag_y[_qp] + 2*0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*_mag_z[_qp] + _alpha[_qp]*_Ms[_qp]*(-2*0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*_mag_x[_qp]*_mag_y[_qp] + 2*0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*_mag_x[_qp]*_mag_z[_qp] + 0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*(Utility::pow<2>(_mag_x[_qp]) - Utility::pow<2>(_mag_y[_qp]) + 3*Utility::pow<2>(_mag_z[_qp])))))*_phi[_j][_qp]*_test[_i][_qp])/
   ((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
    }
    else if (jvar == _disp_x_var)
    {
      return (_g0[_qp]*(4*_B11[_qp]*_grad_phi[_j][_qp](0)*_mag_x[_qp]*(-(_alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_y[_qp]) + _mag_z[_qp]) + _B22[_qp]*(_grad_phi[_j][_qp](2)*(-Utility::pow<2>(_mag_x[_qp]) - 2*_alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_y[_qp]*_mag_z[_qp] + Utility::pow<2>(_mag_z[_qp])) + _grad_phi[_j][_qp](1)*(_mag_y[_qp]*_mag_z[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*(Utility::pow<2>(_mag_x[_qp]) - Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))))*_test[_i][_qp])/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
    }
    else if (jvar == _disp_y_var)
    {
      return (_g0[_qp]*(_B22[_qp]*_mag_y[_qp]*(-(_grad_phi[_j][_qp](2)*_mag_x[_qp]) + _grad_phi[_j][_qp](0)*_mag_z[_qp]) + _alpha[_qp]*_Ms[_qp]*(4*_B11[_qp]*_grad_phi[_j][_qp](1)*_mag_y[_qp]*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_z[_qp])) + _B22[_qp]*(_grad_phi[_j][_qp](0)*_mag_x[_qp] + _grad_phi[_j][_qp](2)*_mag_z[_qp])*(Utility::pow<2>(_mag_x[_qp]) - Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp]))))*_test[_i][_qp])/(2.*(1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
    }
    else if (jvar == _disp_z_var)
    {
      return -(_g0[_qp]*(4*_B11[_qp]*_grad_phi[_j][_qp](2)*_mag_z[_qp]*(_mag_x[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp]) + _B22[_qp]*(_mag_x[_qp]*(_grad_phi[_j][_qp](0)*_mag_x[_qp] + _grad_phi[_j][_qp](1)*_mag_y[_qp]) + _alpha[_qp]*_Ms[_qp]*(2*_grad_phi[_j][_qp](0)*_mag_x[_qp]*_mag_y[_qp] + _grad_phi[_j][_qp](1)*(-Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp])))*_mag_z[_qp] - _grad_phi[_j][_qp](0)*Utility::pow<2>(_mag_z[_qp]) - _alpha[_qp]*_grad_phi[_j][_qp](1)*_Ms[_qp]*Utility::pow<3>(_mag_z[_qp])))*_test[_i][_qp])/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
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
      return (_g0[_qp]*(2*_B11[_qp]*(-_disp_x_grad[_qp](0) + _disp_y_grad[_qp](1))*_mag_y[_qp] + 4*_alpha[_qp]*_B11[_qp]*(-_disp_x_grad[_qp](0) + _disp_z_grad[_qp](2))*_Ms[_qp]*_mag_x[_qp]*_mag_z[_qp] + _B22[_qp]*(0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_z[_qp] + 2*0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*(_mag_x[_qp] - _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp]) + _alpha[_qp]*_Ms[_qp]*(2*0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*_mag_x[_qp]*_mag_y[_qp] + 0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*(3*Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) - Utility::pow<2>(_mag_z[_qp])))))*_phi[_j][_qp]*_test[_i][_qp])/
   ((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
    }
    else if (jvar == _mag_y_var)
    {
      return -((_g0[_qp]*(2*_B11[_qp]*(_disp_x_grad[_qp](0) - _disp_y_grad[_qp](1))*_mag_x[_qp] + 4*_alpha[_qp]*_B11[_qp]*(_disp_y_grad[_qp](1) - _disp_z_grad[_qp](2))*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp] + _B22[_qp]*(2*0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*_mag_y[_qp] + 0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*_mag_z[_qp] + _alpha[_qp]*_Ms[_qp]*(-2*0.5*(_disp_x_grad[_qp](2)+_disp_z_grad[_qp](0))*_mag_x[_qp]*_mag_y[_qp] + 2*0.5*(_disp_x_grad[_qp](1)+_disp_y_grad[_qp](0))*_mag_x[_qp]*_mag_z[_qp] + 0.5*(_disp_y_grad[_qp](2)+_disp_z_grad[_qp](1))*(-Utility::pow<2>(_mag_x[_qp]) - 3*Utility::pow<2>(_mag_y[_qp]) + Utility::pow<2>(_mag_z[_qp])))))*_phi[_j][_qp]*_test[_i][_qp])/
     ((1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp])));
    }
    else if (jvar == _disp_x_var)
    {
      return -(_g0[_qp]*(4*_B11[_qp]*_grad_phi[_j][_qp](0)*_mag_x[_qp]*(_mag_y[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_z[_qp]) + _B22[_qp]*(_grad_phi[_j][_qp](1)*(-Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) + 2*_alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_y[_qp]*_mag_z[_qp]) + _grad_phi[_j][_qp](2)*(_mag_y[_qp]*_mag_z[_qp] - _alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) - Utility::pow<2>(_mag_z[_qp])))))*_test[_i][_qp])/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
    }
    else if (jvar == _disp_y_var)
    {
      return (_g0[_qp]*(4*_B11[_qp]*_grad_phi[_j][_qp](1)*_mag_y[_qp]*(_mag_x[_qp] - _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*_mag_z[_qp]) + _B22[_qp]*(_grad_phi[_j][_qp](0)*(Utility::pow<2>(_mag_x[_qp]) - Utility::pow<2>(_mag_y[_qp]) - 2*_alpha[_qp]*_Ms[_qp]*_mag_x[_qp]*_mag_y[_qp]*_mag_z[_qp]) + _grad_phi[_j][_qp](2)*(_mag_x[_qp]*_mag_z[_qp] + _alpha[_qp]*_Ms[_qp]*_mag_y[_qp]*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) - Utility::pow<2>(_mag_z[_qp])))))*_test[_i][_qp])/
   (2.*(1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
    }
    else if (jvar == _disp_z_var)
    {
      return (_g0[_qp]*(_B22[_qp]*(_grad_phi[_j][_qp](1)*_mag_x[_qp] - _grad_phi[_j][_qp](0)*_mag_y[_qp])*_mag_z[_qp] + _alpha[_qp]*_Ms[_qp]*(4*_B11[_qp]*_grad_phi[_j][_qp](2)*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]))*_mag_z[_qp] + _B22[_qp]*(_grad_phi[_j][_qp](0)*_mag_x[_qp] + _grad_phi[_j][_qp](1)*_mag_y[_qp])*(Utility::pow<2>(_mag_x[_qp]) + Utility::pow<2>(_mag_y[_qp]) - Utility::pow<2>(_mag_z[_qp]))))*_test[_i][_qp])/(2.*(1 + Utility::pow<2>(_alpha[_qp]))*Utility::pow<2>(_Ms[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
