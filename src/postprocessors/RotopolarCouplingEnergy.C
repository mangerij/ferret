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

#include "RotopolarCouplingEnergy.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<RotopolarCouplingEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over the fourth order coupling energy density between AFD and polarization fields.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt vector");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the antiferrodistortive tilt vector");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodistortive tilt vector");
  params.addRequiredParam<Real>("t11", "The coupling constants");
  params.addRequiredParam<Real>("t12", "The coupling constants");
  params.addRequiredParam<Real>("t44", "The coupling constants");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

RotopolarCouplingEnergy::RotopolarCouplingEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
  _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
  _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
  _t11(getParam<Real>("t11")),
  _t12(getParam<Real>("t12")),
  _t44(getParam<Real>("t44")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
RotopolarCouplingEnergy::computeQpIntegral()
{
  return (-((Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t11) - 
   (Utility::pow<2>(_antiferrodis_A_z[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp])) + Utility::pow<2>(_antiferrodis_A_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp])) + Utility::pow<2>(_antiferrodis_A_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])))*_t12 - (_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + _antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t44) * Utility::pow<3>(_len_scale);
}
