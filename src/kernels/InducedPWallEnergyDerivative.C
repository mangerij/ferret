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

#include "InducedPWallEnergyDerivative.h"

registerMooseObject("FerretApp", InducedPWallEnergyDerivative);

template<>
InputParameters validParams<InducedPWallEnergyDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("induced_polar_x", "The x component of the induced polarization");
  params.addRequiredCoupledVar("induced_polar_y", "The y component of the induced polarization");
  params.addCoupledVar("induced_polar_z", 0.0, "The z component of the induced polarization");
  params.addRequiredParam<Real>("G110", "Domain wall coefficient");
  params.addRequiredParam<Real>("G11_G110", "Domain wall coefficient ratio");
  params.addRequiredParam<Real>("G12_G110", "Domain wall coefficient ratio");
  params.addRequiredParam<Real>("G44_G110", "Domain wall coefficient ratio");
  params.addRequiredParam<Real>("G44P_G110", "Domain wall coefficient ratio");
  ///params.set<bool>("use_displaced_mesh") = false;
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

InducedPWallEnergyDerivative::InducedPWallEnergyDerivative(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _induced_polar_x_var(coupled("induced_polar_x")),
  _induced_polar_y_var(coupled("induced_polar_y")),
  _induced_polar_z_var(coupled("induced_polar_z")),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z")),
  _induced_polar_x_grad(coupledGradient("induced_polar_x")),
  _induced_polar_y_grad(coupledGradient("induced_polar_y")),
  _induced_polar_z_grad(coupledGradient("induced_polar_z")),
  _G110(getParam<Real>("G110")),
  _G11(getParam<Real>("G11_G110") * _G110),
  _G12(getParam<Real>("G12_G110") * _G110),
  _G44(getParam<Real>("G44_G110") * _G110),
  _G44P(getParam<Real>("G44P_G110") * _G110),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
InducedPWallEnergyDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return _grad_test[_i][_qp](1)*_G44*(_polar_x_grad[_qp](1) + _polar_y_grad[_qp](0)) + _grad_test[_i][_qp](2)*_G44*(_polar_x_grad[_qp](2) + _polar_z_grad[_qp](0)) + _grad_test[_i][_qp](0)*(_G11*_polar_x_grad[_qp](0) + _G12*(_polar_y_grad[_qp](1) + _polar_z_grad[_qp](2)));
  }
  else if (_component == 1)
  {
    return _grad_test[_i][_qp](0)*_G44*(_polar_x_grad[_qp](1) + _polar_y_grad[_qp](0)) + _grad_test[_i][_qp](2)*_G44*(_polar_y_grad[_qp](2) + _polar_z_grad[_qp](1)) + _grad_test[_i][_qp](1)*(_G11*_polar_y_grad[_qp](1) + _G12*(_polar_x_grad[_qp](0) + _polar_z_grad[_qp](2)));
  }
  else if (_component == 2)
  {
    return _grad_test[_i][_qp](0)*_G44*(_polar_x_grad[_qp](2) + _polar_z_grad[_qp](0)) + _grad_test[_i][_qp](1)*_G44*(_polar_y_grad[_qp](2) + _polar_z_grad[_qp](1)) + _grad_test[_i][_qp](2)*(_G12*(_polar_x_grad[_qp](0) + _polar_y_grad[_qp](1)) + _G11*_polar_z_grad[_qp](2));
  }
  else
    return 0.0;
}

Real
InducedPWallEnergyDerivative::computeQpJacobian()
{
  return 0.0;
}

Real
InducedPWallEnergyDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
