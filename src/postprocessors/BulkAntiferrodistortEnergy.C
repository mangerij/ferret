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

#include "BulkAntiferrodistortEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", BulkAntiferrodistortEnergy);

template<>
InputParameters validParams<BulkAntiferrodistortEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over the local sixth order in the AFD field energy density.");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt vector");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the antiferrodistortive tilt vector");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodistortive tilt vector");
  params.addRequiredParam<Real>("beta1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta11", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta12", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta123", "The coefficients of the Landau expansion");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

BulkAntiferrodistortEnergy::BulkAntiferrodistortEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
  _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
  _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
  _beta1(getParam<Real>("beta1")),
  _beta11(getParam<Real>("beta11")),
  _beta12(getParam<Real>("beta12")),
  _beta111(getParam<Real>("beta111")),
  _beta112(getParam<Real>("beta112")),
  _beta123(getParam<Real>("beta123")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
BulkAntiferrodistortEnergy::computeQpIntegral()
{
  return (
    _beta1 * (Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))
  + _beta11 * (Utility::pow<4>(_antiferrodis_A_x[_qp]) + Utility::pow<4>(_antiferrodis_A_y[_qp]) + Utility::pow<4>(_antiferrodis_A_z[_qp]))
    + _beta12 * (Utility::pow<2>(_antiferrodis_A_x[_qp]) * Utility::pow<2>(_antiferrodis_A_y[_qp])+
	      Utility::pow<2>(_antiferrodis_A_y[_qp]) * Utility::pow<2>(_antiferrodis_A_z[_qp])+
	      Utility::pow<2>(_antiferrodis_A_x[_qp]) * Utility::pow<2>(_antiferrodis_A_z[_qp]))+
    _beta111 * (Utility::pow<6>(_antiferrodis_A_x[_qp]) + Utility::pow<6>(_antiferrodis_A_y[_qp]) + Utility::pow<6>(_antiferrodis_A_z[_qp]))+
    _beta112 * (Utility::pow<4>(_antiferrodis_A_x[_qp]) * (Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))
	      + Utility::pow<4>(_antiferrodis_A_y[_qp]) * (Utility::pow<2>(_antiferrodis_A_z[_qp]) + Utility::pow<2>(_antiferrodis_A_x[_qp]))
	      + Utility::pow<4>(_antiferrodis_A_z[_qp]) * (Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp])))+
	  _beta123 * (pow(_antiferrodis_A_x[_qp], 2) * Utility::pow<2>(_antiferrodis_A_y[_qp]) * Utility::pow<2>(_antiferrodis_A_z[_qp]))) * Utility::pow<3>(_len_scale);
}
