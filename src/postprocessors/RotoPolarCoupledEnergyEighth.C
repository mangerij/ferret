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

#include "RotoPolarCoupledEnergyEighth.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", RotoPolarCoupledEnergyEighth);

InputParameters RotoPolarCoupledEnergyEighth::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over the eighth order coupling energy density between AFD and polarization fields.");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the AFD vector field");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the AFD vector field");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the AFD vector field");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

RotoPolarCoupledEnergyEighth::RotoPolarCoupledEnergyEighth(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
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
   _t1111(getMaterialProperty<Real>("t1111")),
   _t1122(getMaterialProperty<Real>("t1122")),
   _t1212(getMaterialProperty<Real>("t1212")),
   _t42111111(getMaterialProperty<Real>("t42111111")),
   _t24111111(getMaterialProperty<Real>("t24111111")),
   _t42111122(getMaterialProperty<Real>("t42111122")),
   _t24112222(getMaterialProperty<Real>("t24112222")),
   _t42112233(getMaterialProperty<Real>("t42112233")),
   _t24112233(getMaterialProperty<Real>("t24112233")),
   _t42112211(getMaterialProperty<Real>("t42112211")),
   _t24111122(getMaterialProperty<Real>("t24111122")),
   _t42111212(getMaterialProperty<Real>("t42111212")),
   _t42123312(getMaterialProperty<Real>("t42123312")),
   _t24121112(getMaterialProperty<Real>("t24121112")),
   _t24121233(getMaterialProperty<Real>("t24121233")),
   _t6211111111(getMaterialProperty<Real>("t6211111111")),
   _t2611111111(getMaterialProperty<Real>("t2611111111")),
   _t6211111122(getMaterialProperty<Real>("t6211111122")),
   _t2611222222(getMaterialProperty<Real>("t2611222222")),
   _t4411111111(getMaterialProperty<Real>("t4411111111")),
   _t4411112222(getMaterialProperty<Real>("t4411112222")),
   _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
RotoPolarCoupledEnergyEighth::computeQpIntegral()
{
  return _energy_scale*((Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t1111[_qp] + 
   ((Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<2>(_polar_z[_qp]))*_t1122[_qp] + (_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + _antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t1212[_qp] + 
   (Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111111[_qp] + 
   (Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] + 
   ((Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + (Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) + (Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]))*Utility::pow<2>(_polar_z[_qp]))*_t24112222[_qp] + 
   (Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112233[_qp] + 
   (Utility::pow<3>(_antiphase_A_z[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp])*_polar_z[_qp] + Utility::pow<3>(_antiphase_A_y[_qp])*_polar_y[_qp]*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]) + Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp]*(_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t24121112[_qp] + (_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] + 
   (Utility::pow<6>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + Utility::pow<6>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<6>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t2611111111[_qp] + 
   ((Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<6>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + (Utility::pow<6>(_antiphase_A_x[_qp]) + Utility::pow<6>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) + (Utility::pow<6>(_antiphase_A_x[_qp]) + Utility::pow<6>(_antiphase_A_y[_qp]))*Utility::pow<2>(_polar_z[_qp]))*_t2611222222[_qp] + 
   (Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t42111111[_qp] + 
   ((Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<4>(_polar_x[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<4>(_polar_y[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<4>(_polar_z[_qp]))*_t42111122[_qp] + 
   ((_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp])*Utility::pow<3>(_polar_z[_qp]) + Utility::pow<3>(_polar_x[_qp])*(_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_z[_qp]) + Utility::pow<3>(_polar_y[_qp])*(_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_z[_qp]))*_t42111212[_qp] + 
   (Utility::pow<2>(_antiphase_A_z[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])))*_t42112211[_qp] + 
   (Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t42112233[_qp] + 
   (_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] + (Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t4411111111[_qp] + 
   ((Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]))*Utility::pow<4>(_polar_x[_qp]) + (Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]))*Utility::pow<4>(_polar_y[_qp]) + (Utility::pow<4>(_antiphase_A_x[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp]))*Utility::pow<4>(_polar_z[_qp]))*_t4411112222[_qp] + 
   (Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<6>(_polar_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<6>(_polar_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<6>(_polar_z[_qp]))*_t6211111111[_qp] + 
   ((Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<6>(_polar_x[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<6>(_polar_y[_qp]) + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<6>(_polar_z[_qp]))*_t6211111122[_qp]);
}
