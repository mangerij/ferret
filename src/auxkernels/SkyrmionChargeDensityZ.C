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

#include "SkyrmionChargeDensityZ.h"
template<>

InputParameters validParams<SkyrmionChargeDensityZ>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("polar_x_norm", "The x component of the normalized polarization");
  params.addRequiredCoupledVar("polar_y_norm", "The y component of the normalized polarization");
  params.addRequiredCoupledVar("polar_z_norm", "The z component of the normalized polarization");
  return params;
}


SkyrmionChargeDensityZ::SkyrmionChargeDensityZ(const InputParameters & parameters) :
  AuxKernel(parameters),
   _polar_x_norm(coupledValue("polar_x_norm")),
   _polar_y_norm(coupledValue("polar_y_norm")),
   _polar_z_norm(coupledValue("polar_z_norm")),
   _polar_x_norm_grad(coupledGradient("polar_x_norm")),
   _polar_y_norm_grad(coupledGradient("polar_y_norm")),
   _polar_z_norm_grad(coupledGradient("polar_z_norm"))
{
}

Real
SkyrmionChargeDensityZ::computeValue()

{
    return std::abs((1.0 / (4.0 * 3.14159)) * (_polar_z_norm[_qp] * (-_polar_x_norm_grad[_qp](1) * _polar_y_norm_grad[_qp](0) + _polar_x_norm_grad[_qp](0) * _polar_y_norm_grad[_qp](1)) + _polar_y_norm[_qp] * (_polar_x_norm_grad[_qp](1) * _polar_z_norm_grad[_qp](0) - _polar_x_norm_grad[_qp](0) * _polar_z_norm_grad[_qp](1)) + _polar_x_norm[_qp] *(-_polar_y_norm_grad[_qp](1) * _polar_z_norm_grad[_qp](0) + _polar_y_norm_grad[_qp](0) * _polar_z_norm_grad[_qp](1))));
}
