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

#include "RotoBulkEnergyDensity.h"

registerMooseObject("FerretApp", RotoBulkEnergyDensity);

InputParameters RotoBulkEnergyDensity::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the afd vector field");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the afd vector field");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the afd vector field");
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
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
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
  return (_beta123*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + _beta1*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + _beta12*(Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp])) + _beta11*(Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + _beta1123*(Utility::pow<6>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<4>(_antiphase_A_z[_qp])) + _beta1122*(Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<4>(_antiphase_A_z[_qp])) + _beta111*(Utility::pow<6>(_antiphase_A_x[_qp]) + Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<6>(_antiphase_A_z[_qp])) + _beta1111*(Utility::pow<8>(_antiphase_A_x[_qp]) + Utility::pow<8>(_antiphase_A_y[_qp]) + Utility::pow<8>(_antiphase_A_z[_qp])) + 
   _beta112*((Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))) + 
   _beta1112*((Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<6>(_antiphase_A_z[_qp]) + Utility::pow<6>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<6>(_antiphase_A_x[_qp])*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])))) * std::pow(_len_scale,3);
}
