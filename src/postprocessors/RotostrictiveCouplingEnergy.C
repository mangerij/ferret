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

   You should have received a co_polar_y[_qp] of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "RotostrictiveCouplingEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", RotostrictiveCouplingEnergy);

template<>
InputParameters validParams<RotostrictiveCouplingEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over the coupling free energy density between the AFD and elastic fields.");
  params.addRequiredCoupledVar("u_x", "The x component of the displacement");
  params.addRequiredCoupledVar("u_y", "The y component of the displacement");
  params.addCoupledVar("u_z", 0.0, "The z component of the displacement");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt vector");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the antiferrodistortive tilt vector");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodistortive tilt vector");
  params.addRequiredParam<Real>("r11", "The coupling constants");
  params.addRequiredParam<Real>("r12", "The coupling constants");
  params.addRequiredParam<Real>("r44", "The coupling constants");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

RotostrictiveCouplingEnergy::RotostrictiveCouplingEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _u_x_var(coupled("u_x")),
   _u_y_var(coupled("u_y")),
   _u_z_var(coupled("u_z")),
   _antiferrodis_A_x_var(coupled("antiferrodis_A_x")),
   _antiferrodis_A_y_var(coupled("antiferrodis_A_y")),
   _antiferrodis_A_z_var(coupled("antiferrodis_A_z")),
   _u_x_grad(coupledGradient("u_x")),
   _u_y_grad(coupledGradient("u_y")),
   _u_z_grad(coupledGradient("u_z")),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
   _r11(getParam<Real>("r11")),
   _r12(getParam<Real>("r12")),
   _r44(getParam<Real>("r44")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
RotostrictiveCouplingEnergy::computeQpIntegral()
{
  return -(-2.0*_r44*((_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0)))/2.0 + (_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*(_u_x_grad[_qp](2) + _u_z_grad[_qp](0)))/2.0 + (_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)))/2.0) - _r12*((Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))*_u_x_grad[_qp](0) + (Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))*_u_y_grad[_qp](1) + (Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*_u_z_grad[_qp](2)) - 
   _r11*(Utility::pow<2>(_antiferrodis_A_x[_qp])*_u_x_grad[_qp](0) + Utility::pow<2>(_antiferrodis_A_y[_qp])*_u_y_grad[_qp](1) + Utility::pow<2>(_antiferrodis_A_z[_qp])*_u_z_grad[_qp](2))) * Utility::pow<3>(_len_scale);

}
