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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "WindingNumberDensity.h"

registerMooseObject("FerretApp", WindingNumberDensity);

template<>
InputParameters validParams<WindingNumberDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("norm_polar_x", "The x component of the normalized polarization");
  params.addRequiredCoupledVar("norm_polar_y", "The y component of the normalized polarization");
  params.addCoupledVar("norm_polar_z", 0.0, "The z component of the normalized polarization");
  return params;
}

WindingNumberDensity::WindingNumberDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
  _norm_polar_x(coupledValue("norm_polar_x")),
  _norm_polar_y(coupledValue("norm_polar_y")),
  _norm_polar_z(coupledValue("norm_polar_z")),
  _norm_polar_x_grad(coupledGradient("norm_polar_x")),
  _norm_polar_y_grad(coupledGradient("norm_polar_y")),
  _norm_polar_z_grad(coupledGradient("norm_polar_z"))
{
}

Real
WindingNumberDensity::computeValue()
{
  // Pz (-xPy yPx + xPx yPy) + Py (xPz yPx - xPx yPz) +  Px (-xPz yPy + xPy yPz)
  return (1.0 / (4.0 * 3.14159265359)) * (
   _norm_polar_x[_qp] * (-_norm_polar_y_grad[_qp](0) * _norm_polar_x_grad[_qp](1) + _norm_polar_x_grad[_qp](0) * _norm_polar_y_grad[_qp](1))
 + _norm_polar_y[_qp] * (_norm_polar_z_grad[_qp](0) * _norm_polar_x_grad[_qp](1) - _norm_polar_x_grad[_qp](0) * _norm_polar_z_grad[_qp](1))
 + _norm_polar_z[_qp] * (-_norm_polar_z_grad[_qp](0) * _norm_polar_y_grad[_qp](1) + _norm_polar_y_grad[_qp](0) * _norm_polar_z_grad[_qp](1)));
}
