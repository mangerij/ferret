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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "MagHStrongSublatticesCart.h"

class MagHStrongSublatticesCart;

registerMooseObject("FerretApp", MagHStrongSublatticesCart);

InputParameters MagHStrongSublatticesCart::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution for bound magnetic charge (div M1 + div M2)");
  params.addRequiredCoupledVar("mag1_x", "The x component of a constrained sublattice magnetic vector");
  params.addRequiredCoupledVar("mag1_y", "The y component of a constrained sublattice magnetic vector");
  params.addRequiredCoupledVar("mag1_z", "The z component of a constrained sublattice magnetic vector");
  params.addRequiredCoupledVar("mag2_x", "The x component of a constrained sublattice magnetic vector");
  params.addRequiredCoupledVar("mag2_y", "The y component of a constrained sublattice magnetic vector");
  params.addRequiredCoupledVar("mag2_z", "The z component of a constrained sublattice magnetic vector");
  return params;
}

MagHStrongSublatticesCart::MagHStrongSublatticesCart(const InputParameters & parameters)
  :Kernel(parameters),
   _mag1_x_var(coupled("mag1_x")),
   _mag1_y_var(coupled("mag1_y")),
   _mag1_z_var(coupled("mag1_z")),
   _mag1_x(coupledValue("mag1_x")),
   _mag1_y(coupledValue("mag1_y")),
   _mag1_z(coupledValue("mag1_z")),
   _mag1_x_grad(coupledGradient("mag1_x")),
   _mag1_y_grad(coupledGradient("mag1_y")),
   _mag1_z_grad(coupledGradient("mag1_z")),
   _mag2_x_var(coupled("mag2_x")),
   _mag2_y_var(coupled("mag2_y")),
   _mag2_z_var(coupled("mag2_z")),
   _mag2_x(coupledValue("mag2_x")),
   _mag2_y(coupledValue("mag2_y")),
   _mag2_z(coupledValue("mag2_z")),
   _mag2_x_grad(coupledGradient("mag2_x")),
   _mag2_y_grad(coupledGradient("mag2_y")),
   _mag2_z_grad(coupledGradient("mag2_z")),
   _mu0(getMaterialProperty<Real>("mu0")),
   _Ms(getMaterialProperty<Real>("Ms"))
{
}

Real
MagHStrongSublatticesCart::computeQpResidual()
{
  return -0.5*_Ms[_qp]*_mu0[_qp]*(_grad_test[_i][_qp](0)*_mag1_x[_qp]+_grad_test[_i][_qp](1)*_mag1_y[_qp]+_grad_test[_i][_qp](2)*_mag1_z[_qp] + _grad_test[_i][_qp](0)*_mag2_x[_qp]+_grad_test[_i][_qp](1)*_mag2_y[_qp]+_grad_test[_i][_qp](2)*_mag2_z[_qp]);
}

Real
MagHStrongSublatticesCart::computeQpJacobian()
{
  return 0.0;
}

Real
MagHStrongSublatticesCart::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _mag1_x_var)
  {
    return -0.5*_Ms[_qp]*_mu0[_qp]*(_grad_test[_i][_qp](0)*_phi[_j][_qp]);
  }
  else if (jvar == _mag1_y_var)
  {
    return -0.5*_Ms[_qp]*_mu0[_qp]*(_grad_test[_i][_qp](1)*_phi[_j][_qp]);
  }
  else if (jvar == _mag1_z_var)
  {
    return -0.5*_Ms[_qp]*_mu0[_qp]*(_grad_test[_i][_qp](2)*_phi[_j][_qp]);
  }
  else if (jvar == _mag2_x_var)
  {
    return -0.5*_Ms[_qp]*_mu0[_qp]*(_grad_test[_i][_qp](0)*_phi[_j][_qp]);
  }
  else if (jvar == _mag2_y_var)
  {
    return -0.5*_Ms[_qp]*_mu0[_qp]*(_grad_test[_i][_qp](1)*_phi[_j][_qp]);
  }
  else if (jvar == _mag2_z_var)
  {
    return -0.5*_Ms[_qp]*_mu0[_qp]*(_grad_test[_i][_qp](2)*_phi[_j][_qp]);
  }
  else
    return 0.0;
}
