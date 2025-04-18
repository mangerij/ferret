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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "WallEnergyDerivative.h"

registerMooseObject("FerretApp", WallEnergyDerivative);

InputParameters WallEnergyDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to the variation w.r.t polarization of the gradient energy. This Kernel needs to be used in conjunction with Wall2EnergyDerivative!");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  ///params.set<bool>("use_displaced_mesh") = false;
  return params;
}

WallEnergyDerivative::WallEnergyDerivative(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _polar_x_var(coupled("polar_x")),
  _polar_y_var(coupled("polar_y")),
  _polar_z_var(coupled("polar_z")),
  _polar_i_grad((_component==0)? coupledGradient("polar_x") :(_component==1)? coupledGradient("polar_y"): coupledGradient("polar_z")),
  _polar_j_grad((_component==0)? coupledGradient("polar_y"): (_component==1)? coupledGradient("polar_z"): coupledGradient("polar_x")),
  _polar_k_grad((_component==0)? coupledGradient("polar_z"): (_component==1)? coupledGradient("polar_x"): coupledGradient("polar_y")),
  _ii(_component),
  _jj((_component==0)? 1 : (_component==1)? 2: 0),
  _kk((_component==0)? 2 : (_component==1)? 0: 1),
  _G110(getMaterialProperty<Real>("G110")),
  _G11(getMaterialProperty<Real>("G11_G110")),
  _G12(getMaterialProperty<Real>("G12_G110")),
  _G44(getMaterialProperty<Real>("G44_G110")),
  _G44P(getMaterialProperty<Real>("G44P_G110"))
{
}

Real
WallEnergyDerivative::computeQpResidual()
{
  Real Rwall = 0.0;

  Rwall += (_G11[_qp] * _polar_i_grad[_qp](_ii) * _grad_test[_i][_qp](_ii) +
    _G12[_qp] * (_polar_j_grad[_qp](_jj) + _polar_k_grad[_qp](_kk)) * _grad_test[_i][_qp](_ii) +
    _G44[_qp] * (_polar_i_grad[_qp](_jj) + _polar_j_grad[_qp](_ii)) * _grad_test[_i][_qp](_jj) + _G44[_qp] * (_polar_i_grad[_qp](_kk)+_polar_k_grad[_qp](_ii)) * _grad_test[_i][_qp](_kk) +
	   _G44P[_qp] * (_polar_i_grad[_qp](_jj) - _polar_j_grad[_qp](_ii)) * _grad_test[_i][_qp](_jj) + _G44P[_qp] * (_polar_i_grad[_qp](_kk) - _polar_k_grad[_qp](_ii)) * _grad_test[_i][_qp](_kk));
  ///  Moose::out << "\n R_wall-"; std::cout << _component << " = " << Rwall;
  return _G110[_qp]*Rwall;
}

Real
WallEnergyDerivative::computeQpJacobian()
{
  return _G110[_qp]*((_G11[_qp] * _grad_phi[_j][_qp](_ii) * _grad_test[_i][_qp](_ii) +
    (_G44[_qp] + _G44P[_qp]) * _grad_phi[_j][_qp](_jj) * _grad_test[_i][_qp](_jj) +
	  (_G44[_qp] + _G44P[_qp]) * _grad_phi[_j][_qp](_kk) * _grad_test[_i][_qp](_kk)) );
}

Real
WallEnergyDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  //mooseAssert(jvar!=variable().number(),"Something wrong: OffDiag coupled to itself.");
  if(jvar==_polar_x_var || jvar==_polar_y_var || jvar==_polar_z_var)
  {
    const unsigned int _jj = (jvar==_polar_x_var)? 0: (jvar==_polar_y_var)? 1 : 2;
    return _G110[_qp]*((_G12[_qp] * _grad_phi[_j][_qp](_jj) * _grad_test[_i][_qp](_ii) + (_G44[_qp] - _G44P[_qp]) * _grad_phi[_j][_qp](_ii) * _grad_test[_i][_qp](_jj)) );
  }
  else
  {
    return 0.0;
  }
}
