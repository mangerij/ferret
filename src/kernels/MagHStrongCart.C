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

#include "MagHStrongCart.h"

class MagHStrongCart;

registerMooseObject("FerretApp", MagHStrongCart);

template<>
InputParameters validParams<MagHStrongCart>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution for bound magnetic charge (div M)");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetic vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained magnetic vector");
  return params;
}

MagHStrongCart::MagHStrongCart(const InputParameters & parameters)
  :Kernel(parameters),
   _mag_x_var(coupled("mag_x")),
   _mag_y_var(coupled("mag_y")),
   _mag_z_var(coupled("mag_z")),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _mag_x_grad(coupledGradient("mag_x")),
   _mag_y_grad(coupledGradient("mag_y")),
   _mag_z_grad(coupledGradient("mag_z")),
   _mu0(getMaterialProperty<Real>("mu0")),
   _Ms(getMaterialProperty<Real>("Ms"))
{
}

Real
MagHStrongCart::computeQpResidual()
{
  return -_Ms[_qp]*(_grad_test[_i][_qp](0)*_mag_x[_qp]+_grad_test[_i][_qp](1)*_mag_y[_qp]+_grad_test[_i][_qp](2)*_mag_z[_qp]);
}
Real
MagHStrongCart::computeQpJacobian()
{
  return 0.0;
}

Real
MagHStrongCart::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _mag_x_var)
  {
    return -_Ms[_qp]*(_grad_test[_i][_qp](0)*_phi[_j][_qp]);
  }
  else if (jvar == _mag_y_var)
  {
    return -_Ms[_qp]*(_grad_test[_i][_qp](1)*_phi[_j][_qp]);
  }
  else if (jvar == _mag_z_var)
  {
    return -_Ms[_qp]*(_grad_test[_i][_qp](2)*_phi[_j][_qp]);
  }
  else
    return 0.0;
}
