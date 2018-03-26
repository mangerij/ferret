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

#include "RotoBulkEnergyDensity.h"

registerMooseObject("FerretApp", RotoBulkEnergyDensity);

template<>
InputParameters validParams<RotoBulkEnergyDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the afd vector field");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the afd vector field");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the afd vector field");
  params.addRequiredParam<Real>("beta1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta11", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta12", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta123", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta1111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta1112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta1122", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("beta1123", "The coefficients of the Landau expansion");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

RotoBulkEnergyDensity::RotoBulkEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
   _beta1(getParam<Real>("beta1")),
   _beta11(getParam<Real>("beta11")),
   _beta12(getParam<Real>("beta12")),
   _beta111(getParam<Real>("beta111")),
   _beta112(getParam<Real>("beta112")),
   _beta123(getParam<Real>("beta123")),
   _beta1111(getParam<Real>("beta1111")),
   _beta1112(getParam<Real>("beta1112")),
   _beta1122(getParam<Real>("beta1122")),
   _beta1123(getParam<Real>("beta1123")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
RotoBulkEnergyDensity::computeValue()
{
  return (_beta123*std::pow(_antiferrodis_A_x[_qp],2.0)*std::pow(_antiferrodis_A_y[_qp],2.0)*std::pow(_antiferrodis_A_z[_qp],2.0) + _beta1*(std::pow(_antiferrodis_A_x[_qp],2.0) + std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2.0)) + _beta12*(std::pow(_antiferrodis_A_x[_qp],2.0)*std::pow(_antiferrodis_A_y[_qp],2.0) + std::pow(_antiferrodis_A_x[_qp],2.0)*std::pow(_antiferrodis_A_z[_qp],2.0) + std::pow(_antiferrodis_A_y[_qp],2.0)*std::pow(_antiferrodis_A_z[_qp],2.0)) + _beta11*(std::pow(_antiferrodis_A_x[_qp],4) + std::pow(_antiferrodis_A_y[_qp],4.0) + std::pow(_antiferrodis_A_z[_qp],4.0)) + _beta1123*(std::pow(_antiferrodis_A_x[_qp],6)*std::pow(_antiferrodis_A_z[_qp],2.0) + std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],4.0)*std::pow(_antiferrodis_A_z[_qp],2.0) + std::pow(_antiferrodis_A_x[_qp],2.0)*std::pow(_antiferrodis_A_y[_qp],2.0)*std::pow(_antiferrodis_A_z[_qp],4)) + _beta1122*(std::pow(_antiferrodis_A_x[_qp],4)*std::pow(_antiferrodis_A_y[_qp],4.0) + std::pow(_antiferrodis_A_x[_qp],4.0)*std::pow(_antiferrodis_A_z[_qp],4.0) + std::pow(_antiferrodis_A_y[_qp],4)*std::pow(_antiferrodis_A_z[_qp],4.0)) + _beta111*(std::pow(_antiferrodis_A_x[_qp],6.0) + std::pow(_antiferrodis_A_y[_qp],6.0) + std::pow(_antiferrodis_A_z[_qp],6.0)) + _beta1111*(std::pow(_antiferrodis_A_x[_qp],8) + std::pow(_antiferrodis_A_y[_qp],8.0) + std::pow(_antiferrodis_A_z[_qp],8.0)) + 
   _beta112*((std::pow(_antiferrodis_A_x[_qp],2.0) + std::pow(_antiferrodis_A_y[_qp],2))*std::pow(_antiferrodis_A_z[_qp],4) + std::pow(_antiferrodis_A_y[_qp],4.0)*(std::pow(_antiferrodis_A_x[_qp],2.0) + std::pow(_antiferrodis_A_z[_qp],2.0)) + std::pow(_antiferrodis_A_x[_qp],4.0)*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2.0))) + 
   _beta1112*((std::pow(_antiferrodis_A_x[_qp],2.0) + std::pow(_antiferrodis_A_y[_qp],2))*std::pow(_antiferrodis_A_z[_qp],6) + std::pow(_antiferrodis_A_y[_qp],6.0)*(std::pow(_antiferrodis_A_x[_qp],2.0) + std::pow(_antiferrodis_A_z[_qp],2.0)) + std::pow(_antiferrodis_A_x[_qp],6.0)*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2.0)))) * std::pow(_len_scale,3);
}
