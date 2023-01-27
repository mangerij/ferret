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

#include "Wall2EnergyDerivative.h"

registerMooseObject("FerretApp", Wall2EnergyDerivative);

InputParameters Wall2EnergyDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to the variation w.r.t polarization of the gradient energy. This Kernel needs to be used in conjunction with WallEnergyDerivative!");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in OP space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

Wall2EnergyDerivative::Wall2EnergyDerivative(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _polar_x_var(coupled("polar_x")),
  _polar_y_var(coupled("polar_y")),
  _polar_z_var(coupled("polar_z")),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z")),
  _G110(getMaterialProperty<Real>("G110")),
  _G11(getMaterialProperty<Real>("G11_G110")),
  _G12(getMaterialProperty<Real>("G12_G110")),
  _G44(getMaterialProperty<Real>("G44_G110")),
  _G44P(getMaterialProperty<Real>("G44P_G110"))
{
}

Real
Wall2EnergyDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return (_polar_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_G11[_qp] + (_polar_y_grad[_qp](1)*_grad_test[_i][_qp](0) + _polar_z_grad[_qp](2)*_grad_test[_i][_qp](0))*_G12[_qp] + 
   ((2*(_polar_x_grad[_qp](1) + _polar_y_grad[_qp](0))*_grad_test[_i][_qp](1) + 2*(_polar_x_grad[_qp](2) + _polar_z_grad[_qp](0))*_grad_test[_i][_qp](2))*_G44[_qp])/2. + 
   ((2*(_polar_x_grad[_qp](1) - _polar_y_grad[_qp](0))*_grad_test[_i][_qp](1) + 2*(_polar_x_grad[_qp](2) - _polar_z_grad[_qp](0))*_grad_test[_i][_qp](2))*_G44P[_qp])/2.);
  }
  else if (_component == 1)
  {
    return (_polar_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_G11[_qp] + (_polar_x_grad[_qp](0)*_grad_test[_i][_qp](1) + _polar_z_grad[_qp](2)*_grad_test[_i][_qp](1))*_G12[_qp] + 
   ((2*(_polar_x_grad[_qp](1) + _polar_y_grad[_qp](0))*_grad_test[_i][_qp](0) + 2*(_polar_y_grad[_qp](2) + _polar_z_grad[_qp](1))*_grad_test[_i][_qp](2))*_G44[_qp])/2. + 
   ((-2*(_polar_x_grad[_qp](1) - _polar_y_grad[_qp](0))*_grad_test[_i][_qp](0) - 2*(-_polar_y_grad[_qp](2) + _polar_z_grad[_qp](1))*_grad_test[_i][_qp](2))*_G44P[_qp])/2.);
  }
  else if (_component == 2)
  {
    return (_polar_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_G11[_qp] + (_polar_x_grad[_qp](0)*_grad_test[_i][_qp](2) + _polar_y_grad[_qp](1)*_grad_test[_i][_qp](2))*_G12[_qp] + 
   ((2*(_polar_x_grad[_qp](2) + _polar_z_grad[_qp](0))*_grad_test[_i][_qp](0) + 2*(_polar_y_grad[_qp](2) + _polar_z_grad[_qp](1))*_grad_test[_i][_qp](1))*_G44[_qp])/2. + 
   ((-2*(_polar_x_grad[_qp](2) - _polar_z_grad[_qp](0))*_grad_test[_i][_qp](0) + 2*(-_polar_y_grad[_qp](2) + _polar_z_grad[_qp](1))*_grad_test[_i][_qp](1))*_G44P[_qp])/2.);
  }
  else
    return 0.0;
}

Real
Wall2EnergyDerivative::computeQpJacobian()
{
  if (_component == 0)
  {
    return (_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*_G11[_qp] + ((2*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + 2*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*_G44[_qp])/2. + 
   ((2*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + 2*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*_G44P[_qp])/2.);
  }
  else if (_component == 1)
  {
    return (_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*_G11[_qp] + ((2*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + 2*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*_G44[_qp])/2. + 
   ((2*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + 2*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*_G44P[_qp])/2.);
  }
  else if (_component == 2)
  {
    return (_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*_G11[_qp] + ((2*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + 2*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1))*_G44[_qp])/2. + 
   ((2*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + 2*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1))*_G44P[_qp])/2.);
  }
  else
    return 0.0;
}

Real
Wall2EnergyDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return (_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](0)*_G12[_qp] + _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](1)*_G44[_qp] - _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](1)*_G44P[_qp]);
    }
    else if (jvar == _polar_z_var)
    {
      return (_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](0)*_G12[_qp] + _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](2)*_G44[_qp] - _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](2)*_G44P[_qp]);
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
      return (_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](1)*_G12[_qp] + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](0)*_G44[_qp] - _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](0)*_G44P[_qp]);
    }
    else if (jvar == _polar_z_var)
    {
      return (_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](1)*_G12[_qp] + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](2)*_G44[_qp] - _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](2)*_G44P[_qp]);
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
      return (_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](2)*_G12[_qp] + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](0)*_G44[_qp] - _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](0)*_G44P[_qp]);
    }
    else if (jvar == _polar_y_var)
    {
      return (_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](2)*_G12[_qp] + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](1)*_G44[_qp] - _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](1)*_G44P[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
