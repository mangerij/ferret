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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "RotoPolarCouplingEnergyDensity.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", RotoPolarCouplingEnergyDensity);

InputParameters RotoPolarCouplingEnergyDensity::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the afd vector field");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the afd vector field");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the afd vector field");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("t1111", "The coupling coefficients");
  params.addRequiredParam<Real>("t1122", "The coupling coefficients");
  params.addRequiredParam<Real>("t1212", "The coupling coefficients");
  params.addRequiredParam<Real>("t42111111", "The coupling coefficients");
  params.addRequiredParam<Real>("t24111111", "The coupling coefficients");
  params.addRequiredParam<Real>("t42111122", "The coupling coefficients");
  params.addRequiredParam<Real>("t24112222", "The coupling coefficients");
  params.addRequiredParam<Real>("t42112233", "The coupling coefficients");
  params.addRequiredParam<Real>("t24112233", "The coupling coefficients");
  params.addRequiredParam<Real>("t42112211", "The coupling coefficients");
  params.addRequiredParam<Real>("t24111122", "The coupling coefficients");
  params.addRequiredParam<Real>("t42111212", "The coupling coefficients");
  params.addRequiredParam<Real>("t42123312", "The coupling coefficients");
  params.addRequiredParam<Real>("t24121112", "The coupling coefficients");
  params.addRequiredParam<Real>("t24121233", "The coupling coefficients");
  params.addRequiredParam<Real>("t6211111111", "The coupling coefficients");
  params.addRequiredParam<Real>("t2611111111", "The coupling coefficients");
  params.addRequiredParam<Real>("t6211111122", "The coupling coefficients");
  params.addRequiredParam<Real>("t2611222222", "The coupling coefficients");
  params.addRequiredParam<Real>("t4411111111", "The coupling coefficients");
  params.addRequiredParam<Real>("t4411112222", "The coupling coefficients");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

RotoPolarCouplingEnergyDensity::RotoPolarCouplingEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
   _antiphase_A_x_var(coupled("antiphase_A_x")),
   _antiphase_A_y_var(coupled("antiphase_A_y")),
   _antiphase_A_z_var(coupled("antiphase_A_z")),
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _t1111(getParam<Real>("t1111")),
   _t1122(getParam<Real>("t1122")),
   _t1212(getParam<Real>("t1212")),
   _t42111111(getParam<Real>("t42111111")),
   _t24111111(getParam<Real>("t24111111")),
   _t42111122(getParam<Real>("t42111122")),
   _t24112222(getParam<Real>("t24112222")),
   _t42112233(getParam<Real>("t42112233")),
   _t24112233(getParam<Real>("t24112233")),
   _t42112211(getParam<Real>("t42112211")),
   _t24111122(getParam<Real>("t24111122")),
   _t42111212(getParam<Real>("t42111212")),
   _t42123312(getParam<Real>("t42123312")),
   _t24121112(getParam<Real>("t24121112")),
   _t24121233(getParam<Real>("t24121233")),
   _t6211111111(getParam<Real>("t6211111111")),
   _t2611111111(getParam<Real>("t2611111111")),
   _t6211111122(getParam<Real>("t6211111122")),
   _t2611222222(getParam<Real>("t2611222222")),
   _t4411111111(getParam<Real>("t4411111111")),
   _t4411112222(getParam<Real>("t4411112222")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
RotoPolarCouplingEnergyDensity::computeValue()
{
  return (Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t1111 +
   ((Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<2>(_polar_z[_qp]))*_t1122 + (_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + _antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t1212 +
   (Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111111 +
   (Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111122 +
   ((Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + (Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) + (Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]))*Utility::pow<2>(_polar_z[_qp]))*_t24112222 +
   (Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112233 +
   (Utility::pow<3>(_antiphase_A_z[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp])*_polar_z[_qp] + Utility::pow<3>(_antiphase_A_y[_qp])*_polar_y[_qp]*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]) + Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp]*(_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t24121112 + (_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233 +
   (Utility::pow<6>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<6>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<6>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t2611111111 +
   ((Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<6>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + (Utility::pow<6>(_antiphase_A_x[_qp]) + Utility::pow<6>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) + (Utility::pow<6>(_antiphase_A_x[_qp]) + Utility::pow<6>(_antiphase_A_y[_qp]))*Utility::pow<2>(_polar_z[_qp]))*_t2611222222 +
   (Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t42111111 +
   ((Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<4>(_polar_x[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<4>(_polar_y[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<4>(_polar_z[_qp]))*_t42111122 +
   ((_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp])*Utility::pow<3>(_polar_z[_qp]) + Utility::pow<3>(_polar_x[_qp])*(_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_z[_qp]) + Utility::pow<3>(_polar_y[_qp])*(_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_z[_qp]))*_t42111212 +
   (Utility::pow<2>(_antiphase_A_z[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])))*_t42112211 +
   (Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t42112233 +
   (_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312 + (Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t4411111111 +
   ((Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]))*Utility::pow<4>(_polar_x[_qp]) + (Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]))*Utility::pow<4>(_polar_y[_qp]) + (Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]))*Utility::pow<4>(_polar_z[_qp]))*_t4411112222 +
   (Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<6>(_polar_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<6>(_polar_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<6>(_polar_z[_qp]))*_t6211111111 +
   ((Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<6>(_polar_x[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<6>(_polar_y[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<6>(_polar_z[_qp]))*_t6211111122;
}
