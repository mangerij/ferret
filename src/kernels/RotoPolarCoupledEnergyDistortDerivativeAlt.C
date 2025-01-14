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

#include "RotoPolarCoupledEnergyDistortDerivativeAlt.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", RotoPolarCoupledEnergyDistortDerivativeAlt);

InputParameters RotoPolarCoupledEnergyDistortDerivativeAlt::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase tilt");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase tilt");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the antiphase tilt");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

RotoPolarCoupledEnergyDistortDerivativeAlt::RotoPolarCoupledEnergyDistortDerivativeAlt(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
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
   _t4411112222(getMaterialProperty<Real>("t4411112222"))
{
}

Real
RotoPolarCoupledEnergyDistortDerivativeAlt::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (2*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_t1111[_qp] + (2*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 2*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t1122[_qp] +
   (_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp])*_t1212[_qp] + 4*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_t24111111[_qp] +
   (2*_antiphase_A_x[_qp]*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + 2*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) +
      2*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] +
   (4*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 4*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112222[_qp] +
   (2*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112233[_qp] +
   (Utility::pow<3>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_y[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*(_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*
    _t24121112[_qp] + (_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*
    _t24121233[_qp] + 6*Utility::pow<5>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_t2611111111[_qp] +
   (6*Utility::pow<5>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 6*Utility::pow<5>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t2611222222[_qp] +
   2*_antiphase_A_x[_qp]*Utility::pow<4>(_polar_x[_qp])*_t42111111[_qp] + (2*_antiphase_A_x[_qp]*Utility::pow<4>(_polar_y[_qp]) + 2*_antiphase_A_x[_qp]*Utility::pow<4>(_polar_z[_qp]))*_t42111122[_qp] +
   (_antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<3>(_polar_y[_qp]) + _antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<3>(_polar_z[_qp]) + Utility::pow<3>(_polar_x[_qp])*(_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*
    _t42111212[_qp] + 2*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211[_qp] +
   2*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])*_t42112233[_qp] +
   (_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] +
   4*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_x[_qp])*_t4411111111[_qp] +
   (4*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_y[_qp]) + 4*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t4411112222[_qp] +
   2*_antiphase_A_x[_qp]*Utility::pow<6>(_polar_x[_qp])*_t6211111111[_qp] + (2*_antiphase_A_x[_qp]*Utility::pow<6>(_polar_y[_qp]) + 2*_antiphase_A_x[_qp]*Utility::pow<6>(_polar_z[_qp]))*_t6211111122[_qp]);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (2*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp])*_t1111[_qp] + (2*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t1122[_qp] +
   (_antiphase_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t1212[_qp] + 4*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_t24111111[_qp] +
   (2*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2*_antiphase_A_y[_qp]*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) +
      2*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] +
   (4*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 4*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112222[_qp] +
   (2*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24112233[_qp] +
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_y[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_y[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp]*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*
    _t24121112[_qp] + (_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*
    _t24121233[_qp] + 6*Utility::pow<5>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_t2611111111[_qp] +
   (6*Utility::pow<5>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 6*Utility::pow<5>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t2611222222[_qp] +
   2*_antiphase_A_y[_qp]*Utility::pow<4>(_polar_y[_qp])*_t42111111[_qp] + (2*_antiphase_A_y[_qp]*Utility::pow<4>(_polar_x[_qp]) + 2*_antiphase_A_y[_qp]*Utility::pow<4>(_polar_z[_qp]))*_t42111122[_qp] +
   (_antiphase_A_x[_qp]*Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_y[_qp]*Utility::pow<3>(_polar_z[_qp]) + Utility::pow<3>(_polar_y[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*
    _t42111212[_qp] + 2*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211[_qp] +
   2*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp])*_t42112233[_qp] +
   (_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] +
   4*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_y[_qp])*_t4411111111[_qp] +
   (4*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_x[_qp]) + 4*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t4411112222[_qp] +
   2*_antiphase_A_y[_qp]*Utility::pow<6>(_polar_y[_qp])*_t6211111111[_qp] + (2*_antiphase_A_y[_qp]*Utility::pow<6>(_polar_x[_qp]) + 2*_antiphase_A_y[_qp]*Utility::pow<6>(_polar_z[_qp]))*_t6211111122[_qp]);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (2*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp])*_t1111[_qp] + (2*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t1122[_qp] +
   (_antiphase_A_x[_qp]*_polar_x[_qp]*_polar_z[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t1212[_qp] + 4*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp])*_t24111111[_qp] +
   (2*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]) +
      2*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] +
   (4*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 4*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t24112222[_qp] +
   (2*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t24112233[_qp] +
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<3>(_antiphase_A_y[_qp])*_polar_y[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_z[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp])*_polar_z[_qp])*
    _t24121112[_qp] + (2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*
    _t24121233[_qp] + 6*Utility::pow<5>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp])*_t2611111111[_qp] +
   (6*Utility::pow<5>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 6*Utility::pow<5>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t2611222222[_qp] +
   2*_antiphase_A_z[_qp]*Utility::pow<4>(_polar_z[_qp])*_t42111111[_qp] + (2*_antiphase_A_z[_qp]*Utility::pow<4>(_polar_x[_qp]) + 2*_antiphase_A_z[_qp]*Utility::pow<4>(_polar_y[_qp]))*_t42111122[_qp] +
   (_antiphase_A_x[_qp]*Utility::pow<3>(_polar_x[_qp])*_polar_z[_qp] + _antiphase_A_y[_qp]*Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + (_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp])*Utility::pow<3>(_polar_z[_qp]))*
    _t42111212[_qp] + 2*_antiphase_A_z[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp])*_t42112211[_qp] +
   2*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*_t42112233[_qp] +
   (_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp])*_t42123312[_qp] +
   4*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_z[_qp])*_t4411111111[_qp] +
   (4*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_x[_qp]) + 4*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_y[_qp]))*_t4411112222[_qp] +
   2*_antiphase_A_z[_qp]*Utility::pow<6>(_polar_z[_qp])*_t6211111111[_qp] + (2*_antiphase_A_z[_qp]*Utility::pow<6>(_polar_x[_qp]) + 2*_antiphase_A_z[_qp]*Utility::pow<6>(_polar_y[_qp]))*_t6211111122[_qp]);
  }
  else
    return 0.0;
}

Real
RotoPolarCoupledEnergyDistortDerivativeAlt::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*Utility::pow<2>(_polar_x[_qp])*_t1111[_qp] + (2*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp]))*_t1122[_qp] +
   12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_t24111111[_qp] +
   (2*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) +
      2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] +
   (12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112222[_qp] +
   (2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112233[_qp] +
   6*_antiphase_A_x[_qp]*_polar_x[_qp]*(_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp])*_t24121112[_qp] + 2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp]*_t24121233[_qp] +
   30*Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_t2611111111[_qp] +
   (30*Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 30*Utility::pow<4>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t2611222222[_qp] +
   2*Utility::pow<4>(_polar_x[_qp])*_t42111111[_qp] + (2*Utility::pow<4>(_polar_y[_qp]) + 2*Utility::pow<4>(_polar_z[_qp]))*_t42111122[_qp] +
   2*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211[_qp] +
   2*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])*_t42112233[_qp] + 12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_x[_qp])*_t4411111111[_qp] +
   (12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_y[_qp]) + 12*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t4411112222[_qp] +
   2*Utility::pow<6>(_polar_x[_qp])*_t6211111111[_qp] + (2*Utility::pow<6>(_polar_y[_qp]) + 2*Utility::pow<6>(_polar_z[_qp]))*_t6211111122[_qp]);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*Utility::pow<2>(_polar_y[_qp])*_t1111[_qp] + (2*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_z[_qp]))*_t1122[_qp] +
   12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_t24111111[_qp] +
   (2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) +
      2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] +
   (12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112222[_qp] +
   (2*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112233[_qp] +
   6*_antiphase_A_y[_qp]*_polar_y[_qp]*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp])*_t24121112[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp]*_t24121233[_qp] +
   30*Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_t2611111111[_qp] +
   (30*Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 30*Utility::pow<4>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t2611222222[_qp] +
   2*Utility::pow<4>(_polar_y[_qp])*_t42111111[_qp] + (2*Utility::pow<4>(_polar_x[_qp]) + 2*Utility::pow<4>(_polar_z[_qp]))*_t42111122[_qp] +
   2*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211[_qp] +
   2*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp])*_t42112233[_qp] + 12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_y[_qp])*_t4411111111[_qp] +
   (12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_x[_qp]) + 12*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t4411112222[_qp] +
   2*Utility::pow<6>(_polar_y[_qp])*_t6211111111[_qp] + (2*Utility::pow<6>(_polar_x[_qp]) + 2*Utility::pow<6>(_polar_z[_qp]))*_t6211111122[_qp]);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*Utility::pow<2>(_polar_z[_qp])*_t1111[_qp] + (2*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_y[_qp]))*_t1122[_qp] +
   12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp])*_t24111111[_qp] +
   (2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) +
      2*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] +
   (12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t24112222[_qp] +
   (2*Utility::pow<2>(_antiphase_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiphase_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t24112233[_qp] +
   6*_antiphase_A_z[_qp]*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp])*_polar_z[_qp]*_t24121112[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_t24121233[_qp] +
   30*Utility::pow<4>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_z[_qp])*_t2611111111[_qp] +
   (30*Utility::pow<4>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 30*Utility::pow<4>(_antiphase_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t2611222222[_qp] +
   2*Utility::pow<4>(_polar_z[_qp])*_t42111111[_qp] + (2*Utility::pow<4>(_polar_x[_qp]) + 2*Utility::pow<4>(_polar_y[_qp]))*_t42111122[_qp] +
   2*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp])*_t42112211[_qp] +
   2*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*_t42112233[_qp] + 12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_z[_qp])*_t4411111111[_qp] +
   (12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_x[_qp]) + 12*Utility::pow<2>(_antiphase_A_z[_qp])*Utility::pow<4>(_polar_y[_qp]))*_t4411112222[_qp] +
   2*Utility::pow<6>(_polar_z[_qp])*_t6211111111[_qp] + (2*Utility::pow<6>(_polar_x[_qp]) + 2*Utility::pow<6>(_polar_y[_qp]))*_t6211111122[_qp]);
  }
  else
    return 0.0;
}

Real
RotoPolarCoupledEnergyDistortDerivativeAlt::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiphase_A_x[_qp]*_polar_x[_qp]*_t1111[_qp] + (_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp])*_t1212[_qp] + 8*Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp]*_t24111111[_qp] +
   4*_antiphase_A_x[_qp]*(Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_polar_x[_qp]*_t24111122[_qp] +
   (Utility::pow<3>(_antiphase_A_y[_qp])*_polar_y[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_x[_qp])*(_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t24121112[_qp] +
   (_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp] + Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_z[_qp])*_t24121233[_qp] +
   12*Utility::pow<5>(_antiphase_A_x[_qp])*_polar_x[_qp]*_t2611111111[_qp] + 8*_antiphase_A_x[_qp]*Utility::pow<3>(_polar_x[_qp])*_t42111111[_qp] +
   (_antiphase_A_y[_qp]*Utility::pow<3>(_polar_y[_qp]) + _antiphase_A_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + 3*Utility::pow<2>(_polar_x[_qp])*(_antiphase_A_y[_qp]*_polar_y[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t42111212[_qp] +
   4*_antiphase_A_x[_qp]*_polar_x[_qp]*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211[_qp] +
   (_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] +
   16*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_polar_x[_qp])*_t4411111111[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_polar_x[_qp])*_t6211111111[_qp]);
    }
    else if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiphase_A_x[_qp]*_polar_y[_qp]*_t1122[_qp] + _antiphase_A_y[_qp]*_polar_x[_qp]*_t1212[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp]*_t24111122[_qp] +
   8*Utility::pow<3>(_antiphase_A_x[_qp])*_polar_y[_qp]*_t24112222[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp]*_t24112233[_qp] +
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_x[_qp] + Utility::pow<3>(_antiphase_A_y[_qp])*_polar_x[_qp])*_t24121112[_qp] +
   (_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_z[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_x[_qp])*_polar_y[_qp]*_t2611222222[_qp] +
   8*_antiphase_A_x[_qp]*Utility::pow<3>(_polar_y[_qp])*_t42111122[_qp] + (_antiphase_A_y[_qp]*Utility::pow<3>(_polar_x[_qp]) + 3*_antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t42111212[_qp] +
   4*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_t42112211[_qp] + 4*_antiphase_A_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112233[_qp] +
   (2*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] +
   16*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_polar_y[_qp])*_t4411112222[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_polar_y[_qp])*_t6211111122[_qp]);
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiphase_A_x[_qp]*_polar_z[_qp]*_t1122[_qp] + _antiphase_A_z[_qp]*_polar_x[_qp]*_t1212[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_z[_qp]*_t24111122[_qp] +
   8*Utility::pow<3>(_antiphase_A_x[_qp])*_polar_z[_qp]*_t24112222[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_z[_qp]*_t24112233[_qp] +
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_x[_qp])*_t24121112[_qp] +
   (Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_x[_qp])*_polar_z[_qp]*_t2611222222[_qp] +
   8*_antiphase_A_x[_qp]*Utility::pow<3>(_polar_z[_qp])*_t42111122[_qp] + (_antiphase_A_z[_qp]*Utility::pow<3>(_polar_x[_qp]) + 3*_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42111212[_qp] +
   4*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp]*_t42112211[_qp] + 4*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42112233[_qp] +
   (_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 2*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312[_qp] +
   16*Utility::pow<3>(_antiphase_A_x[_qp])*Utility::pow<3>(_polar_z[_qp])*_t4411112222[_qp] + 12*_antiphase_A_x[_qp]*Utility::pow<5>(_polar_z[_qp])*_t6211111122[_qp]);
    }
    else if (jvar == _antiphase_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (_polar_x[_qp]*_polar_y[_qp]*_t1212[_qp] + (4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t24111122[_qp] +
   4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t24112233[_qp] +
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_y[_qp] + 3*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_y[_qp])*_t24121112[_qp] +
   (Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + 2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] +
   (Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + _polar_x[_qp]*Utility::pow<3>(_polar_y[_qp]))*_t42111212[_qp] + _polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42123312[_qp]);
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (_polar_x[_qp]*_polar_z[_qp]*_t1212[_qp] + (4*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] +
   4*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_t24112233[_qp] +
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_z[_qp])*_t24121112[_qp] +
   (2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] +
   (Utility::pow<3>(_polar_x[_qp])*_polar_z[_qp] + _polar_x[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] + _polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42123312[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _polar_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiphase_A_y[_qp]*_polar_x[_qp]*_t1122[_qp] + _antiphase_A_x[_qp]*_polar_y[_qp]*_t1212[_qp] + 4*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_x[_qp]*_t24111122[_qp] +
   8*Utility::pow<3>(_antiphase_A_y[_qp])*_polar_x[_qp]*_t24112222[_qp] + 4*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_t24112233[_qp] +
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_y[_qp] + 3*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp])*_t24121112[_qp] +
   (_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_z[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_y[_qp])*_polar_x[_qp]*_t2611222222[_qp] +
   8*_antiphase_A_y[_qp]*Utility::pow<3>(_polar_x[_qp])*_t42111122[_qp] + (3*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + _antiphase_A_x[_qp]*Utility::pow<3>(_polar_y[_qp]))*_t42111212[_qp] +
   4*_antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_t42112211[_qp] + 4*_antiphase_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112233[_qp] +
   (2*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] +
   16*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<3>(_polar_x[_qp])*_t4411112222[_qp] + 12*_antiphase_A_y[_qp]*Utility::pow<5>(_polar_x[_qp])*_t6211111122[_qp]);
    }
    else if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiphase_A_y[_qp]*_polar_y[_qp]*_t1111[_qp] + (_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp])*_t1212[_qp] + 8*Utility::pow<3>(_antiphase_A_y[_qp])*_polar_y[_qp]*_t24111111[_qp] +
   4*_antiphase_A_y[_qp]*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_polar_y[_qp]*_t24111122[_qp] +
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_y[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t24121112[_qp] +
   (_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_z[_qp])*_t24121233[_qp] +
   12*Utility::pow<5>(_antiphase_A_y[_qp])*_polar_y[_qp]*_t2611111111[_qp] + 8*_antiphase_A_y[_qp]*Utility::pow<3>(_polar_y[_qp])*_t42111111[_qp] +
   (_antiphase_A_x[_qp]*Utility::pow<3>(_polar_x[_qp]) + _antiphase_A_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + 3*Utility::pow<2>(_polar_y[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_z[_qp]*_polar_z[_qp]))*_t42111212[_qp] +
   4*_antiphase_A_y[_qp]*_polar_y[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211[_qp] +
   (_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + _antiphase_A_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312[_qp] +
   16*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<3>(_polar_y[_qp])*_t4411111111[_qp] + 12*_antiphase_A_y[_qp]*Utility::pow<5>(_polar_y[_qp])*_t6211111111[_qp]);
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiphase_A_y[_qp]*_polar_z[_qp]*_t1122[_qp] + _antiphase_A_z[_qp]*_polar_y[_qp]*_t1212[_qp] + 4*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_z[_qp]*_t24111122[_qp] +
   8*Utility::pow<3>(_antiphase_A_y[_qp])*_polar_z[_qp]*_t24112222[_qp] + 4*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_z[_qp]*_t24112233[_qp] +
   (3*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp] + Utility::pow<3>(_antiphase_A_z[_qp])*_polar_y[_qp])*_t24121112[_qp] +
   (2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_y[_qp])*_polar_z[_qp]*_t2611222222[_qp] +
   8*_antiphase_A_y[_qp]*Utility::pow<3>(_polar_z[_qp])*_t42111122[_qp] + (_antiphase_A_z[_qp]*Utility::pow<3>(_polar_y[_qp]) + 3*_antiphase_A_z[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42111212[_qp] +
   4*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42112211[_qp] + 4*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp]*_t42112233[_qp] +
   (_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312[_qp] +
   16*Utility::pow<3>(_antiphase_A_y[_qp])*Utility::pow<3>(_polar_z[_qp])*_t4411112222[_qp] + 12*_antiphase_A_y[_qp]*Utility::pow<5>(_polar_z[_qp])*_t6211111122[_qp]);
    }
    else if (jvar == _antiphase_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (_polar_x[_qp]*_polar_y[_qp]*_t1212[_qp] + (4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t24111122[_qp] +
   4*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t24112233[_qp] +
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_y[_qp] + 3*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_y[_qp])*_t24121112[_qp] +
   (Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + 2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] +
   (Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + _polar_x[_qp]*Utility::pow<3>(_polar_y[_qp]))*_t42111212[_qp] + _polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42123312[_qp]);
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (_polar_y[_qp]*_polar_z[_qp]*_t1212[_qp] + (4*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 4*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] +
   4*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_t24112233[_qp] +
   (3*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t24121112[_qp] +
   (2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] +
   (Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + _polar_y[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] + Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp]*_t42123312[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _polar_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiphase_A_z[_qp]*_polar_x[_qp]*_t1122[_qp] + _antiphase_A_x[_qp]*_polar_z[_qp]*_t1212[_qp] + 4*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp]*_t24111122[_qp] +
   8*Utility::pow<3>(_antiphase_A_z[_qp])*_polar_x[_qp]*_t24112222[_qp] + 4*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_x[_qp]*_t24112233[_qp] +
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_z[_qp] + 3*_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_z[_qp])*_t24121112[_qp] +
   (2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_y[_qp] + _antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_z[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_z[_qp])*_polar_x[_qp]*_t2611222222[_qp] +
   8*_antiphase_A_z[_qp]*Utility::pow<3>(_polar_x[_qp])*_t42111122[_qp] + (3*_antiphase_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + _antiphase_A_x[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] +
   4*_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112211[_qp] + 4*_antiphase_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_t42112233[_qp] +
   (2*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiphase_A_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp])*_t42123312[_qp] +
   16*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<3>(_polar_x[_qp])*_t4411112222[_qp] + 12*_antiphase_A_z[_qp]*Utility::pow<5>(_polar_x[_qp])*_t6211111122[_qp]);
    }
    else if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiphase_A_z[_qp]*_polar_y[_qp]*_t1122[_qp] + _antiphase_A_y[_qp]*_polar_z[_qp]*_t1212[_qp] + 4*Utility::pow<2>(_antiphase_A_y[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp]*_t24111122[_qp] +
   8*Utility::pow<3>(_antiphase_A_z[_qp])*_polar_y[_qp]*_t24112222[_qp] + 4*Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_z[_qp]*_polar_y[_qp]*_t24112233[_qp] +
   (Utility::pow<3>(_antiphase_A_y[_qp])*_polar_z[_qp] + 3*_antiphase_A_y[_qp]*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_z[_qp])*_t24121112[_qp] +
   (2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_z[_qp])*_t24121233[_qp] + 12*Utility::pow<5>(_antiphase_A_z[_qp])*_polar_y[_qp]*_t2611222222[_qp] +
   8*_antiphase_A_z[_qp]*Utility::pow<3>(_polar_y[_qp])*_t42111122[_qp] + (3*_antiphase_A_y[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _antiphase_A_y[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] +
   4*_antiphase_A_z[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112211[_qp] + 4*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_t42112233[_qp] +
   (_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312[_qp] +
   16*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<3>(_polar_y[_qp])*_t4411112222[_qp] + 12*_antiphase_A_z[_qp]*Utility::pow<5>(_polar_y[_qp])*_t6211111122[_qp]);
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiphase_A_z[_qp]*_polar_z[_qp]*_t1111[_qp] + (_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp])*_t1212[_qp] + 8*Utility::pow<3>(_antiphase_A_z[_qp])*_polar_z[_qp]*_t24111111[_qp] +
   4*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*_antiphase_A_z[_qp]*_polar_z[_qp]*_t24111122[_qp] +
   (Utility::pow<3>(_antiphase_A_x[_qp])*_polar_x[_qp] + Utility::pow<3>(_antiphase_A_y[_qp])*_polar_y[_qp] + 3*Utility::pow<2>(_antiphase_A_z[_qp])*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp]))*_t24121112[_qp] +
   (_antiphase_A_x[_qp]*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_antiphase_A_y[_qp]*_polar_y[_qp])*_t24121233[_qp] +
   12*Utility::pow<5>(_antiphase_A_z[_qp])*_polar_z[_qp]*_t2611111111[_qp] + 8*_antiphase_A_z[_qp]*Utility::pow<3>(_polar_z[_qp])*_t42111111[_qp] +
   (_antiphase_A_x[_qp]*Utility::pow<3>(_polar_x[_qp]) + _antiphase_A_y[_qp]*Utility::pow<3>(_polar_y[_qp]) + 3*(_antiphase_A_x[_qp]*_polar_x[_qp] + _antiphase_A_y[_qp]*_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t42111212[_qp] +
   4*_antiphase_A_z[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*_polar_z[_qp]*_t42112211[_qp] +
   (_antiphase_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + _antiphase_A_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t42123312[_qp] +
   16*Utility::pow<3>(_antiphase_A_z[_qp])*Utility::pow<3>(_polar_z[_qp])*_t4411111111[_qp] + 12*_antiphase_A_z[_qp]*Utility::pow<5>(_polar_z[_qp])*_t6211111111[_qp]);
    }
    else if (jvar == _antiphase_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (_polar_x[_qp]*_polar_z[_qp]*_t1212[_qp] + (4*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] +
   4*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_t24112233[_qp] +
   (3*Utility::pow<2>(_antiphase_A_x[_qp])*_polar_x[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_x[_qp]*_polar_z[_qp])*_t24121112[_qp] +
   (2*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + Utility::pow<2>(_antiphase_A_y[_qp])*_polar_x[_qp]*_polar_z[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] +
   (Utility::pow<3>(_polar_x[_qp])*_polar_z[_qp] + _polar_x[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] + _polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42123312[_qp]);
    }
    else if (jvar == _antiphase_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (_polar_y[_qp]*_polar_z[_qp]*_t1212[_qp] + (4*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 4*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122[_qp] +
   4*_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_t24112233[_qp] +
   (3*Utility::pow<2>(_antiphase_A_y[_qp])*_polar_y[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiphase_A_z[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t24121112[_qp] +
   (2*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + 2*_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<2>(_antiphase_A_x[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t24121233[_qp] +
   (Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + _polar_y[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212[_qp] + Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp]*_t42123312[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
