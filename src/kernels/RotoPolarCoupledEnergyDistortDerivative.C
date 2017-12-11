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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "RotoPolarCoupledEnergyDistortDerivative.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<RotoPolarCoupledEnergyDistortDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the antiferrodistortive tilt");
  params.addRequiredCoupledVar("antiferrodis_A_z", "The z component of the antiferrodistortive tilt");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
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
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

RotoPolarCoupledEnergyDistortDerivative::RotoPolarCoupledEnergyDistortDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _antiferrodis_A_x_var(coupled("antiferrodis_A_x")),
   _antiferrodis_A_y_var(coupled("antiferrodis_A_y")),
   _antiferrodis_A_z_var(coupled("antiferrodis_A_z")),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
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
RotoPolarCoupledEnergyDistortDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_t1111 + (2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t1122 + (_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp])*_t1212 + 4*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_t24111111 + 
   (2*_antiferrodis_A_x[_qp]*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111122 + (4*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 4*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112222 + 
   (2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112233 + (Utility::pow<3>(_antiferrodis_A_y[_qp])*_polar_x[_qp]*_polar_y[_qp] + Utility::pow<3>(_antiferrodis_A_z[_qp])*_polar_x[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiferrodis_A_x[_qp])*_polar_x[_qp]*(_antiferrodis_A_y[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp]))*_t24121112 + 
   (_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233 + 6*Utility::pow<5>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_t2611111111 + (6*Utility::pow<5>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 6*Utility::pow<5>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t2611222222 + 
   2*_antiferrodis_A_x[_qp]*Utility::pow<4>(_polar_x[_qp])*_t42111111 + (2*_antiferrodis_A_x[_qp]*Utility::pow<4>(_polar_y[_qp]) + 2*_antiferrodis_A_x[_qp]*Utility::pow<4>(_polar_z[_qp]))*_t42111122 + (_antiferrodis_A_y[_qp]*_polar_x[_qp]*Utility::pow<3>(_polar_y[_qp]) + _antiferrodis_A_z[_qp]*_polar_x[_qp]*Utility::pow<3>(_polar_z[_qp]) + Utility::pow<3>(_polar_x[_qp])*(_antiferrodis_A_y[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp]))*_t42111212 + 
   2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211 + 2*_antiferrodis_A_x[_qp]*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])*_t42112233 + (_antiferrodis_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312 + 4*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_polar_x[_qp])*_t4411111111 + 
   (4*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 4*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t4411112222 + 2*_antiferrodis_A_x[_qp]*Utility::pow<6>(_polar_x[_qp])*_t6211111111 + (2*_antiferrodis_A_x[_qp]*Utility::pow<6>(_polar_y[_qp]) + 2.0*_antiferrodis_A_x[_qp]*Utility::pow<6>(_polar_z[_qp]))*_t6211111122);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (2*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_y[_qp])*_t1111 + (2.0*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t1122 + (_antiferrodis_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t1212 + 4*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_t24111111 + 
   (2*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2.0*_antiferrodis_A_y[_qp]*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) + 2*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111122 + (4*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 4*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112222 + 
   (2*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2.0*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24112233 + (Utility::pow<3>(_antiferrodis_A_x[_qp])*_polar_x[_qp]*_polar_y[_qp] + Utility::pow<3>(_antiferrodis_A_z[_qp])*_polar_y[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_y[_qp]*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp]))*_t24121112 + 
   (_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + 2.0*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233 + 6.0*Utility::pow<5>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_t2611111111 + (6*Utility::pow<5>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 6*Utility::pow<5>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t2611222222 + 
   2*_antiferrodis_A_y[_qp]*Utility::pow<4>(_polar_y[_qp])*_t42111111 + (2*_antiferrodis_A_y[_qp]*Utility::pow<4>(_polar_x[_qp]) + 2*_antiferrodis_A_y[_qp]*Utility::pow<4>(_polar_z[_qp]))*_t42111122 + (_antiferrodis_A_x[_qp]*Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_y[_qp]*Utility::pow<3>(_polar_z[_qp]) + Utility::pow<3>(_polar_y[_qp])*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp]))*_t42111212 + 
   2*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211 + 2*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp])*_t42112233 + (_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp] + _antiferrodis_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312 + 4*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_polar_y[_qp])*_t4411111111 + 
   (4*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_polar_x[_qp]) + 4*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t4411112222 + 2*_antiferrodis_A_y[_qp]*Utility::pow<6>(_polar_y[_qp])*_t6211111111 + (2*_antiferrodis_A_y[_qp]*Utility::pow<6>(_polar_x[_qp]) + 2*_antiferrodis_A_y[_qp]*Utility::pow<6>(_polar_z[_qp]))*_t6211111122);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (2*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_z[_qp])*_t1111 + (2*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t1122 + (_antiferrodis_A_x[_qp]*_polar_x[_qp]*_polar_z[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t1212 + 4*Utility::pow<3>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_z[_qp])*_t24111111 + 
   (2*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 2*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122 + (4*Utility::pow<3>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 4*Utility::pow<3>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t24112222 + 
   (2*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t24112233 + (Utility::pow<3>(_antiferrodis_A_x[_qp])*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<3>(_antiferrodis_A_y[_qp])*_polar_y[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiferrodis_A_z[_qp])*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp])*_polar_z[_qp])*_t24121112 + 
   (2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233 + 6*Utility::pow<5>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_z[_qp])*_t2611111111 + (6*Utility::pow<5>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 6*Utility::pow<5>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t2611222222 + 
   2*_antiferrodis_A_z[_qp]*Utility::pow<4>(_polar_z[_qp])*_t42111111 + (2*_antiferrodis_A_z[_qp]*Utility::pow<4>(_polar_x[_qp]) + 2*_antiferrodis_A_z[_qp]*Utility::pow<4>(_polar_y[_qp]))*_t42111122 + (_antiferrodis_A_x[_qp]*Utility::pow<3>(_polar_x[_qp])*_polar_z[_qp] + _antiferrodis_A_y[_qp]*Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + (_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp])*Utility::pow<3>(_polar_z[_qp]))*_t42111212 + 
   2*_antiferrodis_A_z[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp])*_t42112211 + 2*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*_t42112233 + (_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp] + _antiferrodis_A_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp])*_t42123312 + 4*Utility::pow<3>(_antiferrodis_A_z[_qp])*Utility::pow<4>(_polar_z[_qp])*_t4411111111 + 
   (4*Utility::pow<3>(_antiferrodis_A_z[_qp])*Utility::pow<4>(_polar_x[_qp]) + 4*Utility::pow<3>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t4411112222 + 2*_antiferrodis_A_z[_qp]*Utility::pow<6>(_polar_z[_qp])*_t6211111111 + (2*_antiferrodis_A_z[_qp]*Utility::pow<6>(_polar_x[_qp]) + 2*_antiferrodis_A_z[_qp]*Utility::pow<6>(_polar_y[_qp]))*_t6211111122);
  }
  else
    return 0.0;
}

Real
RotoPolarCoupledEnergyDistortDerivative::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*Utility::pow<2>(_polar_x[_qp])*_t1111 + (2*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_polar_z[_qp]))*_t1122 + 12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_t24111111 + 
   (2*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111122 + (12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112222 + 
   (2*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112233 + 6*_antiferrodis_A_x[_qp]*_polar_x[_qp]*(_antiferrodis_A_y[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp])*_t24121112 + 2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp]*_t24121233 + 30*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_x[_qp])*_t2611111111 + 
   (30*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 30*Utility::pow<4>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t2611222222 + 2*Utility::pow<4>(_polar_x[_qp])*_t42111111 + (2*Utility::pow<4>(_polar_y[_qp]) + 2*Utility::pow<4>(_polar_z[_qp]))*_t42111122 + 2*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211 + 
   2*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])*_t42112233 + 12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_polar_x[_qp])*_t4411111111 + (12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + 12*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t4411112222 + 2*Utility::pow<6>(_polar_x[_qp])*_t6211111111 + 
   (2*Utility::pow<6>(_polar_y[_qp]) + 2*Utility::pow<6>(_polar_z[_qp]))*_t6211111122);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*Utility::pow<2>(_polar_y[_qp])*_t1111 + (2*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_z[_qp]))*_t1122 + 12*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_t24111111 + 
   (2*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))*Utility::pow<2>(_polar_y[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24111122 + (12*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 12*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112222 + 
   (2*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t24112233 + 6*_antiferrodis_A_y[_qp]*_polar_y[_qp]*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp])*_t24121112 + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp]*_t24121233 + 30*Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_y[_qp])*_t2611111111 + 
   (30*Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 30*Utility::pow<4>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t2611222222 + 2*Utility::pow<4>(_polar_y[_qp])*_t42111111 + (2*Utility::pow<4>(_polar_x[_qp]) + 2*Utility::pow<4>(_polar_z[_qp]))*_t42111122 + 2*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211 + 
   2*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp])*_t42112233 + 12*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_polar_y[_qp])*_t4411111111 + (12*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_polar_x[_qp]) + 12*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<4>(_polar_z[_qp]))*_t4411112222 + 2*Utility::pow<6>(_polar_y[_qp])*_t6211111111 + 
   (2*Utility::pow<6>(_polar_x[_qp]) + 2*Utility::pow<6>(_polar_z[_qp]))*_t6211111122);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*Utility::pow<2>(_polar_z[_qp])*_t1111 + (2*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_polar_y[_qp]))*_t1122 + 12*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_z[_qp])*_t24111111 + 
   (2*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_y[_qp]) + 2*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*Utility::pow<2>(_polar_z[_qp]))*_t24111122 + (12*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 12*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t24112222 + 
   (2*Utility::pow<2>(_antiferrodis_A_y[_qp])*Utility::pow<2>(_polar_x[_qp]) + 2*Utility::pow<2>(_antiferrodis_A_x[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t24112233 + 6*_antiferrodis_A_z[_qp]*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp])*_polar_z[_qp]*_t24121112 + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_t24121233 + 30*Utility::pow<4>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_z[_qp])*_t2611111111 + 
   (30*Utility::pow<4>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_x[_qp]) + 30*Utility::pow<4>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t2611222222 + 2*Utility::pow<4>(_polar_z[_qp])*_t42111111 + (2*Utility::pow<4>(_polar_x[_qp]) + 2*Utility::pow<4>(_polar_y[_qp]))*_t42111122 + 2*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp])*_t42112211 + 
   2*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*_t42112233 + 12*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<4>(_polar_z[_qp])*_t4411111111 + (12*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<4>(_polar_x[_qp]) + 12*Utility::pow<2>(_antiferrodis_A_z[_qp])*Utility::pow<2>(_polar_y[_qp]))*_t4411112222 + 2*Utility::pow<6>(_polar_z[_qp])*_t6211111111 + 
   (2*Utility::pow<6>(_polar_x[_qp]) + 2*Utility::pow<6>(_polar_y[_qp]))*_t6211111122);
  }
  else
    return 0.0;
}

Real
RotoPolarCoupledEnergyDistortDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_x[_qp]*_polar_x[_qp]*_t1111 + (_antiferrodis_A_y[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp])*_t1212 + 8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_polar_x[_qp]*_t24111111 + 4*_antiferrodis_A_x[_qp]*(Utility::pow<2>(_antiferrodis_A_y[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))*_polar_x[_qp]*_t24111122 + (Utility::pow<3>(_antiferrodis_A_y[_qp])*_polar_y[_qp] + Utility::pow<3>(_antiferrodis_A_z[_qp])*_polar_z[_qp] + 3*Utility::pow<2>(_antiferrodis_A_x[_qp])*(_antiferrodis_A_y[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp]))*_t24121112 + 
   (_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_y[_qp] + Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp]*_polar_z[_qp])*_t24121233 + 12*Utility::pow<5>(_antiferrodis_A_x[_qp])*_polar_x[_qp]*_t2611111111 + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_polar_x[_qp])*_t42111111 + (_antiferrodis_A_y[_qp]*Utility::pow<3>(_polar_y[_qp]) + _antiferrodis_A_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + 3*Utility::pow<2>(_polar_x[_qp])*(_antiferrodis_A_y[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp]))*_t42111212 + 
   4*_antiferrodis_A_x[_qp]*_polar_x[_qp]*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211 + (_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312 + 16*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_polar_x[_qp])*_t4411111111 + 12*_antiferrodis_A_x[_qp]*Utility::pow<5>(_polar_x[_qp])*_t6211111111);
    }
    else if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_x[_qp]*_polar_y[_qp]*_t1122 + _antiferrodis_A_y[_qp]*_polar_x[_qp]*_t1212 + 4*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_y[_qp]*_t24111122 + 8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_polar_y[_qp]*_t24112222 + 4*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_y[_qp]*_t24112233 + (3*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_polar_x[_qp] + Utility::pow<3>(_antiferrodis_A_y[_qp])*_polar_x[_qp])*_t24121112 + 
   (_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_x[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_z[_qp])*_t24121233 + 12*Utility::pow<5>(_antiferrodis_A_x[_qp])*_polar_y[_qp]*_t2611222222 + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_polar_y[_qp])*_t42111122 + (_antiferrodis_A_y[_qp]*Utility::pow<3>(_polar_x[_qp]) + 3*_antiferrodis_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t42111212 + 4*_antiferrodis_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_t42112211 + 
   4*_antiferrodis_A_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112233 + (2*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiferrodis_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312 + 8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_polar_y[_qp]*_t4411112222 + 12*_antiferrodis_A_x[_qp]*Utility::pow<5>(_polar_y[_qp])*_t6211111122);
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_x[_qp]*_polar_z[_qp]*_t1122 + _antiferrodis_A_z[_qp]*_polar_x[_qp]*_t1212 + 4*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_z[_qp]*_t24111122 + 8*Utility::pow<3>(_antiferrodis_A_x[_qp])*_polar_z[_qp]*_t24112222 + 4*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_z[_qp]*_t24112233 + (3*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp]*_polar_x[_qp] + Utility::pow<3>(_antiferrodis_A_z[_qp])*_polar_x[_qp])*_t24121112 + 
   (Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp]*_polar_x[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp])*_t24121233 + 12*Utility::pow<5>(_antiferrodis_A_x[_qp])*_polar_z[_qp]*_t2611222222 + 8*_antiferrodis_A_x[_qp]*Utility::pow<3>(_polar_z[_qp])*_t42111122 + (_antiferrodis_A_z[_qp]*Utility::pow<3>(_polar_x[_qp]) + 3*_antiferrodis_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42111212 + 4*_antiferrodis_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp]*_t42112211 + 
   4*_antiferrodis_A_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42112233 + (_antiferrodis_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]) + 2*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312 + 16*Utility::pow<3>(_antiferrodis_A_x[_qp])*Utility::pow<3>(_polar_z[_qp])*_t4411112222 + 12*_antiferrodis_A_x[_qp]*Utility::pow<5>(_polar_z[_qp])*_t6211111122);
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (_polar_x[_qp]*_polar_y[_qp]*_t1212 + (4*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t24111122 + 4*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t24112233 + (3*Utility::pow<2>(_antiferrodis_A_x[_qp])*_polar_x[_qp]*_polar_y[_qp] + 3*Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_x[_qp]*_polar_y[_qp])*_t24121112 + 
   (Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + 2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233 + (Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + _polar_x[_qp]*Utility::pow<3>(_polar_y[_qp]))*_t42111212 + _polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42123312);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (_polar_x[_qp]*_polar_z[_qp]*_t1212 + (4*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122 + 4*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_t24112233 + (3*Utility::pow<2>(_antiferrodis_A_x[_qp])*_polar_x[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_x[_qp]*_polar_z[_qp])*_t24121112 + 
   (2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_x[_qp]*_polar_z[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233 + (Utility::pow<3>(_polar_x[_qp])*_polar_z[_qp] + _polar_x[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212 + _polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42123312);
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
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_t1122 + _antiferrodis_A_x[_qp]*_polar_y[_qp]*_t1212 + 4*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_t24111122 + 8*Utility::pow<3>(_antiferrodis_A_y[_qp])*_polar_x[_qp]*_t24112222 + 4*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_x[_qp]*_t24112233 + (Utility::pow<3>(_antiferrodis_A_x[_qp])*_polar_y[_qp] + 3*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_y[_qp])*_t24121112 + 
   (_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_y[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_z[_qp])*_t24121233 + 12*Utility::pow<5>(_antiferrodis_A_y[_qp])*_polar_x[_qp]*_t2611222222 + 8*_antiferrodis_A_y[_qp]*Utility::pow<3>(_polar_x[_qp])*_t42111122 + (3*_antiferrodis_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + _antiferrodis_A_x[_qp]*Utility::pow<3>(_polar_y[_qp]))*_t42111212 + 4*_antiferrodis_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_t42112211 + 
   4*_antiferrodis_A_y[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112233 + (2*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiferrodis_A_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312 + 16*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_polar_x[_qp])*_t4411112222 + 12*_antiferrodis_A_y[_qp]*Utility::pow<5>(_polar_x[_qp])*_t6211111122);
    }
    else if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_y[_qp]*_polar_y[_qp]*_t1111 + (_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp])*_t1212 + 8*Utility::pow<3>(_antiferrodis_A_y[_qp])*_polar_y[_qp]*_t24111111 + 4*_antiferrodis_A_y[_qp]*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_z[_qp]))*_polar_y[_qp]*_t24111122 + (Utility::pow<3>(_antiferrodis_A_x[_qp])*_polar_x[_qp] + Utility::pow<3>(_antiferrodis_A_z[_qp])*_polar_z[_qp] + 3*Utility::pow<2>(_antiferrodis_A_y[_qp])*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp]))*_t24121112 + 
   (_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_x[_qp] + Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp]*_polar_z[_qp])*_t24121233 + 12*Utility::pow<5>(_antiferrodis_A_y[_qp])*_polar_y[_qp]*_t2611111111 + 8*_antiferrodis_A_y[_qp]*Utility::pow<3>(_polar_y[_qp])*_t42111111 + (_antiferrodis_A_x[_qp]*Utility::pow<3>(_polar_x[_qp]) + _antiferrodis_A_z[_qp]*Utility::pow<3>(_polar_z[_qp]) + 3*Utility::pow<2>(_polar_y[_qp])*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp]))*_t42111212 + 
   4*_antiferrodis_A_y[_qp]*_polar_y[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_t42112211 + (_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + _antiferrodis_A_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42123312 + 16*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_polar_y[_qp])*_t4411111111 + 12*_antiferrodis_A_y[_qp]*Utility::pow<5>(_polar_y[_qp])*_t6211111111);
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_y[_qp]*_polar_z[_qp]*_t1122 + _antiferrodis_A_z[_qp]*_polar_y[_qp]*_t1212 + 4*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_z[_qp]*_t24111122 + 8*Utility::pow<3>(_antiferrodis_A_y[_qp])*_polar_z[_qp]*_t24112222 + 4*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_polar_z[_qp]*_t24112233 + (3*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp]*_polar_y[_qp] + Utility::pow<3>(_antiferrodis_A_z[_qp])*_polar_y[_qp])*_t24121112 + 
   (2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp] + Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp]*_polar_y[_qp])*_t24121233 + 12*Utility::pow<5>(_antiferrodis_A_y[_qp])*_polar_z[_qp]*_t2611222222 + 8*_antiferrodis_A_y[_qp]*Utility::pow<3>(_polar_z[_qp])*_t42111122 + (_antiferrodis_A_z[_qp]*Utility::pow<3>(_polar_y[_qp]) + 3*_antiferrodis_A_z[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t42111212 + 4*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42112211 + 
   4*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp]*_t42112233 + (_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + 2*_antiferrodis_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312 + 16*Utility::pow<3>(_antiferrodis_A_y[_qp])*Utility::pow<3>(_polar_z[_qp])*_t4411112222 + 12*_antiferrodis_A_y[_qp]*Utility::pow<5>(_polar_z[_qp])*_t6211111122);
    }
    else if (jvar == _antiferrodis_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (_polar_x[_qp]*_polar_y[_qp]*_t1212 + (4*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t24111122 + 4*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t24112233 + (3*Utility::pow<2>(_antiferrodis_A_x[_qp])*_polar_x[_qp]*_polar_y[_qp] + 3*Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_x[_qp]*_polar_y[_qp])*_t24121112 + 
   (Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_x[_qp]*_polar_y[_qp] + 2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233 + (Utility::pow<3>(_polar_x[_qp])*_polar_y[_qp] + _polar_x[_qp]*Utility::pow<3>(_polar_y[_qp]))*_t42111212 + _polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42123312);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (_polar_y[_qp]*_polar_z[_qp]*_t1212 + (4*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 4*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122 + 4*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_t24112233 + (3*Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_y[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t24121112 + 
   (2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<2>(_antiferrodis_A_x[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t24121233 + (Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + _polar_y[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212 + Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp]*_t42123312);
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
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_t1122 + _antiferrodis_A_x[_qp]*_polar_z[_qp]*_t1212 + 4*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_t24111122 + 8*Utility::pow<3>(_antiferrodis_A_z[_qp])*_polar_x[_qp]*_t24112222 + 4*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_t24112233 + (Utility::pow<3>(_antiferrodis_A_x[_qp])*_polar_z[_qp] + 3*_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_z[_qp])*_t24121112 + 
   (2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp] + _antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_z[_qp])*_t24121233 + 12*Utility::pow<5>(_antiferrodis_A_z[_qp])*_polar_x[_qp]*_t2611222222 + 8*_antiferrodis_A_z[_qp]*Utility::pow<3>(_polar_x[_qp])*_t42111122 + (3*_antiferrodis_A_x[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + _antiferrodis_A_x[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212 + 4*_antiferrodis_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112211 + 
   4*_antiferrodis_A_z[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_t42112233 + (2*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiferrodis_A_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp])*_t42123312 + 16*Utility::pow<3>(_antiferrodis_A_z[_qp])*Utility::pow<3>(_polar_x[_qp])*_t4411112222 + 12*_antiferrodis_A_z[_qp]*Utility::pow<5>(_polar_x[_qp])*_t6211111122);
    }
    else if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_t1122 + _antiferrodis_A_y[_qp]*_polar_z[_qp]*_t1212 + 4*Utility::pow<2>(_antiferrodis_A_y[_qp])*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_t24111122 + 8*Utility::pow<3>(_antiferrodis_A_z[_qp])*_polar_y[_qp]*_t24112222 + 4*Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_t24112233 + (Utility::pow<3>(_antiferrodis_A_y[_qp])*_polar_z[_qp] + 3*_antiferrodis_A_y[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_z[_qp])*_t24121112 + 
   (2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp] + Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_polar_z[_qp])*_t24121233 + 12*Utility::pow<5>(_antiferrodis_A_z[_qp])*_polar_y[_qp]*_t2611222222 + 8*_antiferrodis_A_z[_qp]*Utility::pow<3>(_polar_y[_qp])*_t42111122 + (3*_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp] + _antiferrodis_A_y[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212 + 4*_antiferrodis_A_z[_qp]*_polar_y[_qp]*Utility::pow<2>(_polar_z[_qp])*_t42112211 + 
   4*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_t42112233 + (_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_z[_qp] + 2*_antiferrodis_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312 + 8*Utility::pow<3>(_antiferrodis_A_z[_qp])*_polar_y[_qp]*_t4411112222 + 12*_antiferrodis_A_z[_qp]*Utility::pow<5>(_polar_y[_qp])*_t6211111122);
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_z[_qp]*_polar_z[_qp]*_t1111 + (_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp])*_t1212 + 8*Utility::pow<3>(_antiferrodis_A_z[_qp])*_polar_z[_qp]*_t24111111 + 4*(Utility::pow<2>(_antiferrodis_A_x[_qp]) + Utility::pow<2>(_antiferrodis_A_y[_qp]))*_antiferrodis_A_z[_qp]*_polar_z[_qp]*_t24111122 + (Utility::pow<3>(_antiferrodis_A_x[_qp])*_polar_x[_qp] + Utility::pow<3>(_antiferrodis_A_y[_qp])*_polar_y[_qp] + 3*Utility::pow<2>(_antiferrodis_A_z[_qp])*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp]))*_t24121112 + 
   (_antiferrodis_A_x[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_x[_qp] + Utility::pow<2>(_antiferrodis_A_x[_qp])*_antiferrodis_A_y[_qp]*_polar_y[_qp])*_t24121233 + 12*Utility::pow<5>(_antiferrodis_A_z[_qp])*_polar_z[_qp]*_t2611111111 + 8*_antiferrodis_A_z[_qp]*Utility::pow<3>(_polar_z[_qp])*_t42111111 + (_antiferrodis_A_x[_qp]*Utility::pow<3>(_polar_x[_qp]) + _antiferrodis_A_y[_qp]*Utility::pow<3>(_polar_y[_qp]) + 3*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]))*_t42111212 + 
   4*_antiferrodis_A_z[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*_polar_z[_qp]*_t42112211 + (_antiferrodis_A_y[_qp]*Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp] + _antiferrodis_A_x[_qp]*_polar_x[_qp]*Utility::pow<2>(_polar_y[_qp]))*_t42123312 + 16*Utility::pow<3>(_antiferrodis_A_z[_qp])*Utility::pow<3>(_polar_z[_qp])*_t4411111111 + 12*_antiferrodis_A_z[_qp]*Utility::pow<5>(_polar_z[_qp])*_t6211111111);
    }
    else if (jvar == _antiferrodis_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (_polar_x[_qp]*_polar_z[_qp]*_t1212 + (4*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_x[_qp]) + 4*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122 + 4*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_y[_qp])*_t24112233 + (3*Utility::pow<2>(_antiferrodis_A_x[_qp])*_polar_x[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_x[_qp]*_polar_z[_qp])*_t24121112 + 
   (2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_x[_qp]*_polar_z[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t24121233 + (Utility::pow<3>(_polar_x[_qp])*_polar_z[_qp] + _polar_x[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212 + _polar_x[_qp]*Utility::pow<2>(_polar_y[_qp])*_polar_z[_qp]*_t42123312);
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (_polar_y[_qp]*_polar_z[_qp]*_t1212 + (4*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_y[_qp]) + 4*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_z[_qp]))*_t24111122 + 4*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*Utility::pow<2>(_polar_x[_qp])*_t24112233 + (3*Utility::pow<2>(_antiferrodis_A_y[_qp])*_polar_y[_qp]*_polar_z[_qp] + 3*Utility::pow<2>(_antiferrodis_A_z[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t24121112 + 
   (2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_z[_qp] + Utility::pow<2>(_antiferrodis_A_x[_qp])*_polar_y[_qp]*_polar_z[_qp])*_t24121233 + (Utility::pow<3>(_polar_y[_qp])*_polar_z[_qp] + _polar_y[_qp]*Utility::pow<3>(_polar_z[_qp]))*_t42111212 + Utility::pow<2>(_polar_x[_qp])*_polar_y[_qp]*_polar_z[_qp]*_t42123312);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
