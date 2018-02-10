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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "MagHStrong.h"

class MagHStrong;

template<>
InputParameters validParams<MagHStrong>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetization");
  params.addRequiredCoupledVar("mag_z", "The z component of the magnetization");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

MagHStrong::MagHStrong(const InputParameters & parameters)
  :Kernel(parameters),
   _mag_x_var(coupled("mag_x")),
   _mag_y_var(coupled("mag_y")),
   _mag_z_var(coupled("mag_z")),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
MagHStrong::computeQpResidual()
{
  return - (_mag_x[_qp] * _grad_test[_i][_qp](0) + _mag_y[_qp] * _grad_test[_i][_qp](1) + _mag_z[_qp] * _grad_test[_i][_qp](2)) * std::pow(_len_scale, 2.0);
}
Real
MagHStrong::computeQpJacobian()
{
  return 0.0;
}

Real
MagHStrong::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _mag_x_var)
    return - _phi[_j][_qp] * _grad_test[_i][_qp](0) * std::pow(_len_scale, 2.0);
  else if (jvar == _mag_y_var)
    return - _phi[_j][_qp] * _grad_test[_i][_qp](1) * std::pow(_len_scale, 2.0);
  else if (jvar == _mag_z_var)
    return - _phi[_j][_qp] * _grad_test[_i][_qp](2) * std::pow(_len_scale, 2.0);
  else
    return 0.0;
}
