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

#include "AFDAntiphaseEnergyDerivative.h"

registerMooseObject("FerretApp", AFDAntiphaseEnergyDerivative);

template<>
InputParameters validParams<AFDAntiphaseEnergyDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates the residual for the local gradients in the antiferrodistortive vector field");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt vector");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the antiferrodistortive tilt vector");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodistortive tilt vector");
  params.addRequiredParam<Real>("H110", "Gradient coefficient in the AFD ordering");
  params.addRequiredParam<Real>("H11_H110", "Gradient coefficient in the AFD ordering");
  params.addRequiredParam<Real>("H12_H110", "Gradient coefficient in the AFD ordering");
  params.addRequiredParam<Real>("H44_H110", "Gradient coefficient in the AFD ordering");
  params.addRequiredParam<Real>("H44P_H110", "Gradient coefficient in the AFD ordering");
  ///params.set<bool>("use_displaced_mesh") = false;
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

AFDAntiphaseEnergyDerivative::AFDAntiphaseEnergyDerivative(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _antiferrodis_A_x_var(coupled("antiferrodis_A_x")),
  _antiferrodis_A_y_var(coupled("antiferrodis_A_y")),
  _antiferrodis_A_z_var(coupled("antiferrodis_A_z")),
  _antiferrodis_A_i_grad((_component==0)? coupledGradient("antiferrodis_A_x") :(_component==1)? coupledGradient("antiferrodis_A_y"): coupledGradient("antiferrodis_A_z")),
  _antiferrodis_A_j_grad((_component==0)? coupledGradient("antiferrodis_A_y"): (_component==1)? coupledGradient("antiferrodis_A_z"): coupledGradient("antiferrodis_A_x")),
  _antiferrodis_A_k_grad((_component==0)? coupledGradient("antiferrodis_A_z"): (_component==1)? coupledGradient("antiferrodis_A_x"): coupledGradient("antiferrodis_A_y")),
  _ii(_component),
  _jj((_component==0)? 1 : (_component==1)? 2: 0),
  _kk((_component==0)? 2 : (_component==1)? 0: 1),
  _H110(getParam<Real>("H110")),
  _H11(getParam<Real>("H11_H110") * _H110),
  _H12(getParam<Real>("H12_H110") * _H110),
  _H44(getParam<Real>("H44_H110") * _H110),
  _H44P(getParam<Real>("H44P_H110") * _H110),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
AFDAntiphaseEnergyDerivative::computeQpResidual()
{
  Real Rwall = 0.0;

  Rwall += (_H11 * _antiferrodis_A_i_grad[_qp](_ii) * _grad_test[_i][_qp](_ii) + _H12 * (_antiferrodis_A_j_grad[_qp](_jj) + _antiferrodis_A_k_grad[_qp](_kk)) * _grad_test[_i][_qp](_ii) +
    _H44 * (_antiferrodis_A_i_grad[_qp](_jj) + _antiferrodis_A_j_grad[_qp](_ii)) * _grad_test[_i][_qp](_jj) + _H44 * (_antiferrodis_A_i_grad[_qp](_kk)+_antiferrodis_A_k_grad[_qp](_ii)) * _grad_test[_i][_qp](_kk) + _H44P * (_antiferrodis_A_i_grad[_qp](_jj) - _antiferrodis_A_j_grad[_qp](_ii)) * _grad_test[_i][_qp](_jj) + _H44P * (_antiferrodis_A_i_grad[_qp](_kk) - _antiferrodis_A_k_grad[_qp](_ii)) * _grad_test[_i][_qp](_kk)) * _len_scale;
  ///  Moose::out << "\n R_wall-"; std::cout << _component << " = " << Rwall;
  return Rwall;
}

Real
AFDAntiphaseEnergyDerivative::computeQpJacobian()
{
  return (_H11 * _grad_phi[_j][_qp](_ii) * _grad_test[_i][_qp](_ii) + (_H44 + _H44P) * _grad_phi[_j][_qp](_jj) * _grad_test[_i][_qp](_jj) + (_H44 + _H44P) * _grad_phi[_j][_qp](_kk) * _grad_test[_i][_qp](_kk)) * _len_scale;
}

Real
AFDAntiphaseEnergyDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if(jvar==_antiferrodis_A_x_var || jvar==_antiferrodis_A_y_var || jvar==_antiferrodis_A_z_var)
  {
    const unsigned int _jj = (jvar==_antiferrodis_A_x_var)? 0: (jvar==_antiferrodis_A_y_var)? 1 : 2;
    return (_H12 * _grad_phi[_j][_qp](_jj) * _grad_test[_i][_qp](_ii) + (_H44 - _H44P) * _grad_phi[_j][_qp](_ii) * _grad_test[_i][_qp](_jj)) * _len_scale;
  }
  else
  {
    return 0.0;
  }
}
