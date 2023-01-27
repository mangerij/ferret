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

#include "AFDWall2EnergyDerivative.h"

registerMooseObject("FerretApp", AFDWall2EnergyDerivative);

InputParameters AFDWall2EnergyDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase_Aization");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase_Aization");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the antiphase_Aization");
  return params;
}

AFDWall2EnergyDerivative::AFDWall2EnergyDerivative(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _antiphase_A_x_var(coupled("antiphase_A_x")),
  _antiphase_A_y_var(coupled("antiphase_A_y")),
  _antiphase_A_z_var(coupled("antiphase_A_z")),
  _antiphase_A_x_grad(coupledGradient("antiphase_A_x")),
  _antiphase_A_y_grad(coupledGradient("antiphase_A_y")),
  _antiphase_A_z_grad(coupledGradient("antiphase_A_z")),
  _H110(getMaterialProperty<Real>("H110")),
  _H11(getMaterialProperty<Real>("H11_H110")),
  _H12(getMaterialProperty<Real>("H12_H110")),
  _H44(getMaterialProperty<Real>("H44_H110")),
  _H44P(getMaterialProperty<Real>("H44P_H110"))
{
}

Real
AFDWall2EnergyDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return (_antiphase_A_x_grad[_qp](0)*_grad_test[_i][_qp](0)*_H11[_qp] + (_antiphase_A_y_grad[_qp](1)*_grad_test[_i][_qp](0) + _antiphase_A_z_grad[_qp](2)*_grad_test[_i][_qp](0))*_H12[_qp] + 
   ((2*(_antiphase_A_x_grad[_qp](1) + _antiphase_A_y_grad[_qp](0))*_grad_test[_i][_qp](1) + 2*(_antiphase_A_x_grad[_qp](2) + _antiphase_A_z_grad[_qp](0))*_grad_test[_i][_qp](2))*_H44[_qp])/2. + 
   ((2*(_antiphase_A_x_grad[_qp](1) - _antiphase_A_y_grad[_qp](0))*_grad_test[_i][_qp](1) + 2*(_antiphase_A_x_grad[_qp](2) - _antiphase_A_z_grad[_qp](0))*_grad_test[_i][_qp](2))*_H44P[_qp])/2.);
  }
  else if (_component == 1)
  {
    return (_antiphase_A_y_grad[_qp](1)*_grad_test[_i][_qp](1)*_H11[_qp] + (_antiphase_A_x_grad[_qp](0)*_grad_test[_i][_qp](1) + _antiphase_A_z_grad[_qp](2)*_grad_test[_i][_qp](1))*_H12[_qp] + 
   ((2*(_antiphase_A_x_grad[_qp](1) + _antiphase_A_y_grad[_qp](0))*_grad_test[_i][_qp](0) + 2*(_antiphase_A_y_grad[_qp](2) + _antiphase_A_z_grad[_qp](1))*_grad_test[_i][_qp](2))*_H44[_qp])/2. + 
   ((-2*(_antiphase_A_x_grad[_qp](1) - _antiphase_A_y_grad[_qp](0))*_grad_test[_i][_qp](0) - 2*(-_antiphase_A_y_grad[_qp](2) + _antiphase_A_z_grad[_qp](1))*_grad_test[_i][_qp](2))*_H44P[_qp])/2.);
  }
  else if (_component == 2)
  {
    return (_antiphase_A_z_grad[_qp](2)*_grad_test[_i][_qp](2)*_H11[_qp] + (_antiphase_A_x_grad[_qp](0)*_grad_test[_i][_qp](2) + _antiphase_A_y_grad[_qp](1)*_grad_test[_i][_qp](2))*_H12[_qp] + 
   ((2*(_antiphase_A_x_grad[_qp](2) + _antiphase_A_z_grad[_qp](0))*_grad_test[_i][_qp](0) + 2*(_antiphase_A_y_grad[_qp](2) + _antiphase_A_z_grad[_qp](1))*_grad_test[_i][_qp](1))*_H44[_qp])/2. + 
   ((-2*(_antiphase_A_x_grad[_qp](2) - _antiphase_A_z_grad[_qp](0))*_grad_test[_i][_qp](0) + 2*(-_antiphase_A_y_grad[_qp](2) + _antiphase_A_z_grad[_qp](1))*_grad_test[_i][_qp](1))*_H44P[_qp])/2.);
  }
  else
    return 0.0;
}

Real
AFDWall2EnergyDerivative::computeQpJacobian()
{
  if (_component == 0)
  {
    return (_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0)*_H11[_qp] + ((2*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + 2*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*_H44[_qp])/2. + 
   ((2*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1) + 2*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*_H44P[_qp])/2.);
  }
  else if (_component == 1)
  {
    return (_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1)*_H11[_qp] + ((2*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + 2*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*_H44[_qp])/2. + 
   ((2*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + 2*_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2))*_H44P[_qp])/2.);
  }
  else if (_component == 2)
  {
    return (_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](2)*_H11[_qp] + ((2*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + 2*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1))*_H44[_qp])/2. + 
   ((2*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0) + 2*_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](1))*_H44P[_qp])/2.);
  }
  else
    return 0.0;
}

Real
AFDWall2EnergyDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiphase_A_y_var)
    {
      return (_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](0)*_H12[_qp] + _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](1)*_H44[_qp] - _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](1)*_H44P[_qp]);
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return (_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](0)*_H12[_qp] + _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](2)*_H44[_qp] - _grad_phi[_j][_qp](0)*_grad_test[_i][_qp](2)*_H44P[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _antiphase_A_x_var)
    {
      return (_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](1)*_H12[_qp] + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](0)*_H44[_qp] - _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](0)*_H44P[_qp]);
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return (_grad_phi[_j][_qp](2)*_grad_test[_i][_qp](1)*_H12[_qp] + _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](2)*_H44[_qp] - _grad_phi[_j][_qp](1)*_grad_test[_i][_qp](2)*_H44P[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _antiphase_A_x_var)
    {
      return (_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](2)*_H12[_qp] + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](0)*_H44[_qp] - _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](0)*_H44P[_qp]);
    }
    else if (jvar == _antiphase_A_y_var)
    {
      return (_grad_phi[_j][_qp](1)*_grad_test[_i][_qp](2)*_H12[_qp] + _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](1)*_H44[_qp] - _grad_phi[_j][_qp](2)*_grad_test[_i][_qp](1)*_H44P[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
