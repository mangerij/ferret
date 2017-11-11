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

#include "RotoPolarCoupledEnergyPolarDerivative.h"
#include<cmath>

template<>
InputParameters validParams<RotoPolarCoupledEnergyPolarDerivative>()
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

RotoPolarCoupledEnergyPolarDerivative::RotoPolarCoupledEnergyPolarDerivative(const InputParameters & parameters)
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
RotoPolarCoupledEnergyPolarDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (2*std::pow(_antiferrodis_A_x[_qp],2)*_polar_x[_qp]*_t1111 + 2*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_polar_x[_qp]*_t1122 + ((_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_y[_qp] + _antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_z[_qp])*_t1212)/4. + 2*std::pow(_antiferrodis_A_x[_qp],4)*_polar_x[_qp]*_t24111111 + (std::pow(_antiferrodis_A_x[_qp],2)*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_polar_x[_qp]*_t24111122)/3. + 2*(std::pow(_antiferrodis_A_y[_qp],4) + std::pow(_antiferrodis_A_z[_qp],4))*_polar_x[_qp]*_t24112222 + (std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_antiferrodis_A_z[_qp],2)*_polar_x[_qp]*_t24112233)/3. + ((_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],3)*_polar_y[_qp] + _antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],3)*_polar_z[_qp] + std::pow(_antiferrodis_A_x[_qp],3)*(_antiferrodis_A_y[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp]))*_t24121112)/8. + 
   ((_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_y[_qp] + _antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*_polar_z[_qp])*_t24121233)/24. + 2*(std::pow(_antiferrodis_A_y[_qp],6) + std::pow(_antiferrodis_A_z[_qp],6))*_polar_x[_qp]*_t2611222222 + 4*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_x[_qp],3)*_t42111111 + 4*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*std::pow(_polar_x[_qp],3)*_t42111122 + ((std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*std::pow(_polar_y[_qp],3) + std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_z[_qp]*std::pow(_polar_z[_qp],3) + 3*std::pow(_polar_x[_qp],2)*(_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_polar_y[_qp] + _antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_z[_qp]))*_t42111212)/8. + ((2*std::pow(_antiferrodis_A_y[_qp],2)*_polar_x[_qp]*std::pow(_polar_y[_qp],2) + 2*std::pow(_antiferrodis_A_z[_qp],2)*_polar_x[_qp]*std::pow(_polar_z[_qp],2) + 2*std::pow(_antiferrodis_A_x[_qp],2)*_polar_x[_qp]*(std::pow(_polar_y[_qp],2) + std::pow(_polar_z[_qp],2)))*_t42112211)/6. + ((2*std::pow(_antiferrodis_A_z[_qp],2)*_polar_x[_qp]*std::pow(_polar_y[_qp],2) + 2*std::pow(_antiferrodis_A_y[_qp],2)*_polar_x[_qp]*std::pow(_polar_z[_qp],2))*_t42112233)/6. + ((2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_y[_qp],2)*_polar_z[_qp] + _antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_y[_qp]*std::pow(_polar_z[_qp],2))*_t42123312)/24. + 4*std::pow(_antiferrodis_A_x[_qp],4)*std::pow(_polar_x[_qp],3)*_t4411111111 + 4*(std::pow(_antiferrodis_A_y[_qp],4) + std::pow(_antiferrodis_A_z[_qp],4))*std::pow(_polar_x[_qp],3)*_t4411112222 + 
   6*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_x[_qp],5)*_t6211111111 + 6*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*std::pow(_polar_x[_qp],5)*_t6211111122);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (2*std::pow(_antiferrodis_A_y[_qp],2)*_polar_y[_qp]*_t1111 + 2*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_polar_y[_qp]*_t1122 + ((_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_z[_qp])*_t1212)/4. + 2*std::pow(_antiferrodis_A_y[_qp],4)*_polar_y[_qp]*_t24111111 + (std::pow(_antiferrodis_A_y[_qp],2)*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_polar_y[_qp]*_t24111122)/3. + 
   2*(std::pow(_antiferrodis_A_x[_qp],4) + std::pow(_antiferrodis_A_z[_qp],4))*_polar_y[_qp]*_t24112222 + (std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_z[_qp],2)*_polar_y[_qp]*_t24112233)/3. + ((std::pow(_antiferrodis_A_x[_qp],3)*_antiferrodis_A_y[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],3)*_polar_z[_qp] + std::pow(_antiferrodis_A_y[_qp],3)*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp]))*_t24121112)/8. + 
   ((_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_x[_qp] + std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_z[_qp])*_t24121233)/24. + 2*(std::pow(_antiferrodis_A_x[_qp],6) + std::pow(_antiferrodis_A_z[_qp],6))*_polar_y[_qp]*_t2611222222 + 4*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_y[_qp],3)*_t42111111 + 4*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*std::pow(_polar_y[_qp],3)*_t42111122 + 
   ((_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_x[_qp],3) + std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*std::pow(_polar_z[_qp],3) + 3*std::pow(_polar_y[_qp],2)*(std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_z[_qp]))*_t42111212)/8. + 
   ((2*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_x[_qp],2)*_polar_y[_qp] + 2*std::pow(_antiferrodis_A_z[_qp],2)*_polar_y[_qp]*std::pow(_polar_z[_qp],2) + 2*std::pow(_antiferrodis_A_y[_qp],2)*_polar_y[_qp]*(std::pow(_polar_x[_qp],2) + std::pow(_polar_z[_qp],2)))*_t42112211)/6. + ((2*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_x[_qp],2)*_polar_y[_qp] + 2*std::pow(_antiferrodis_A_x[_qp],2)*_polar_y[_qp]*std::pow(_polar_z[_qp],2))*_t42112233)/6. + 
   ((_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_x[_qp],2)*_polar_z[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_x[_qp]*std::pow(_polar_z[_qp],2))*_t42123312)/24. + 4*std::pow(_antiferrodis_A_y[_qp],4)*std::pow(_polar_y[_qp],3)*_t4411111111 + 2*(std::pow(_antiferrodis_A_x[_qp],4) + std::pow(_antiferrodis_A_z[_qp],4))*_polar_y[_qp]*_t4411112222 + 
   6*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_y[_qp],5)*_t6211111111 + 6*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*std::pow(_polar_y[_qp],5)*_t6211111122);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (2*std::pow(_antiferrodis_A_z[_qp],2)*_polar_z[_qp]*_t1111 + 2*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*_polar_z[_qp]*_t1122 + ((_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp])*_t1212)/4. + 2*std::pow(_antiferrodis_A_z[_qp],4)*_polar_z[_qp]*_t24111111 + ((std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*std::pow(_antiferrodis_A_z[_qp],2)*_polar_z[_qp]*_t24111122)/3. + 2*(std::pow(_antiferrodis_A_x[_qp],4) + std::pow(_antiferrodis_A_y[_qp],4))*_polar_z[_qp]*_t24112222 + (std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],2)*_polar_z[_qp]*_t24112233)/3. + ((std::pow(_antiferrodis_A_x[_qp],3)*_antiferrodis_A_z[_qp]*_polar_x[_qp] + std::pow(_antiferrodis_A_y[_qp],3)*_antiferrodis_A_z[_qp]*_polar_y[_qp] + std::pow(_antiferrodis_A_z[_qp],3)*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp]))*_t24121112)/8. + 
   ((_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*_polar_x[_qp] + std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp])*_t24121233)/24. + 2*(std::pow(_antiferrodis_A_x[_qp],6) + std::pow(_antiferrodis_A_y[_qp],6))*_polar_z[_qp]*_t2611222222 + 4*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_z[_qp],3)*_t42111111 + 4*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*std::pow(_polar_z[_qp],3)*_t42111122 + ((_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_x[_qp],3) + _antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_y[_qp],3) + 3*(std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_z[_qp]*_polar_x[_qp] + std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*_polar_y[_qp])*std::pow(_polar_z[_qp],2))*_t42111212)/8. + 
   ((2*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_x[_qp],2)*_polar_z[_qp] + 2*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_y[_qp],2)*_polar_z[_qp] + 2*std::pow(_antiferrodis_A_z[_qp],2)*(std::pow(_polar_x[_qp],2) + std::pow(_polar_y[_qp],2))*_polar_z[_qp])*_t42112211)/6. + ((2*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_x[_qp],2)*_polar_z[_qp] + 2*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_y[_qp],2)*_polar_z[_qp])*_t42112233)/6. + ((_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_x[_qp],2)*_polar_y[_qp] + _antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*std::pow(_polar_y[_qp],2) + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312)/24. + 4*std::pow(_antiferrodis_A_z[_qp],4)*std::pow(_polar_z[_qp],3)*_t4411111111 + 4*(std::pow(_antiferrodis_A_x[_qp],4) + std::pow(_antiferrodis_A_y[_qp],4))*std::pow(_polar_z[_qp],3)*_t4411112222 + 
   6*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_z[_qp],5)*_t6211111111 + 6*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*std::pow(_polar_z[_qp],5)*_t6211111122);
  }
  else
    return 0.0;
}

Real
RotoPolarCoupledEnergyPolarDerivative::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*std::pow(_antiferrodis_A_x[_qp],2)*_t1111 + 2*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_t1122 + 2*std::pow(_antiferrodis_A_x[_qp],4)*_t24111111 + (std::pow(_antiferrodis_A_x[_qp],2)*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_t24111122)/3. + 2*(std::pow(_antiferrodis_A_y[_qp],4) + std::pow(_antiferrodis_A_z[_qp],4))*_t24112222 + (std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_antiferrodis_A_z[_qp],2)*_t24112233)/3. + 2*(std::pow(_antiferrodis_A_y[_qp],6) + std::pow(_antiferrodis_A_z[_qp],6))*_t2611222222 + 12*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_x[_qp],2)*_t42111111 + 12*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*std::pow(_polar_x[_qp],2)*_t42111122 + 
   (3*_polar_x[_qp]*(_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_polar_y[_qp] + _antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_z[_qp])*_t42111212)/4. + ((2*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_y[_qp],2) + 2*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_z[_qp],2) + 2*std::pow(_antiferrodis_A_x[_qp],2)*(std::pow(_polar_y[_qp],2) + std::pow(_polar_z[_qp],2)))*_t42112211)/6. + 
   ((2*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_y[_qp],2) + 2*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_z[_qp],2))*_t42112233)/6. + (_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp]*_t42123312)/12. + 12*std::pow(_antiferrodis_A_x[_qp],4)*std::pow(_polar_x[_qp],2)*_t4411111111 + 12*(std::pow(_antiferrodis_A_y[_qp],4) + std::pow(_antiferrodis_A_z[_qp],4))*std::pow(_polar_x[_qp],2)*_t4411112222 + 
   30*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_x[_qp],4)*_t6211111111 + 30*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*std::pow(_polar_x[_qp],4)*_t6211111122);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*std::pow(_antiferrodis_A_y[_qp],2)*_t1111 + 2*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_t1122 + 2*std::pow(_antiferrodis_A_y[_qp],4)*_t24111111 + (std::pow(_antiferrodis_A_y[_qp],2)*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_t24111122)/3. + 2*(std::pow(_antiferrodis_A_x[_qp],4) + std::pow(_antiferrodis_A_z[_qp],4))*_t24112222 + (std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_z[_qp],2)*_t24112233)/3. + 2*(std::pow(_antiferrodis_A_x[_qp],6) + std::pow(_antiferrodis_A_z[_qp],6))*_t2611222222 + 12*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_y[_qp],2)*_t42111111 + 12*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*std::pow(_polar_y[_qp],2)*_t42111122 + 
   (3*_polar_y[_qp]*(std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_z[_qp])*_t42111212)/4. + ((2*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_x[_qp],2) + 2*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_z[_qp],2) + 2*std::pow(_antiferrodis_A_y[_qp],2)*(std::pow(_polar_x[_qp],2) + std::pow(_polar_z[_qp],2)))*_t42112211)/6. + 
   ((2*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_x[_qp],2) + 2*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_z[_qp],2))*_t42112233)/6. + (_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp]*_t42123312)/12. + 12*std::pow(_antiferrodis_A_y[_qp],4)*std::pow(_polar_y[_qp],2)*_t4411111111 + 2*(std::pow(_antiferrodis_A_x[_qp],4) + std::pow(_antiferrodis_A_z[_qp],4))*_t4411112222 + 
   30*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_y[_qp],4)*_t6211111111 + 30*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*std::pow(_polar_y[_qp],4)*_t6211111122);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (2*std::pow(_antiferrodis_A_z[_qp],2)*_t1111 + 2*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*_t1122 + 2*std::pow(_antiferrodis_A_z[_qp],4)*_t24111111 + ((std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*std::pow(_antiferrodis_A_z[_qp],2)*_t24111122)/3. + 2*(std::pow(_antiferrodis_A_x[_qp],4) + std::pow(_antiferrodis_A_y[_qp],4))*_t24112222 + (std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_antiferrodis_A_y[_qp],2)*_t24112233)/3. + 2*(std::pow(_antiferrodis_A_x[_qp],6) + std::pow(_antiferrodis_A_y[_qp],6))*_t2611222222 + 12*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_z[_qp],2)*_t42111111 + 12*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*std::pow(_polar_z[_qp],2)*_t42111122 + 
   (3*(std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_z[_qp]*_polar_x[_qp] + std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*_polar_y[_qp])*_polar_z[_qp]*_t42111212)/4. + ((2*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_x[_qp],2) + 2*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_y[_qp],2) + 2*std::pow(_antiferrodis_A_z[_qp],2)*(std::pow(_polar_x[_qp],2) + std::pow(_polar_y[_qp],2)))*_t42112211)/6. + 
   ((2*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_x[_qp],2) + 2*std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_y[_qp],2))*_t42112233)/6. + (_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_t42123312)/12. + 12*std::pow(_antiferrodis_A_z[_qp],4)*std::pow(_polar_z[_qp],2)*_t4411111111 + 12*(std::pow(_antiferrodis_A_x[_qp],4) + std::pow(_antiferrodis_A_y[_qp],4))*std::pow(_polar_z[_qp],2)*_t4411112222 + 
   30*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_z[_qp],4)*_t6211111111 + 30*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*std::pow(_polar_z[_qp],4)*_t6211111122);
  }
  else
    return 0.0;
}

Real
RotoPolarCoupledEnergyPolarDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  ((_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_t1212)/4. + ((std::pow(_antiferrodis_A_x[_qp],3)*_antiferrodis_A_y[_qp] + _antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],3))*_t24121112)/8. + (_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_t24121233)/24. + ((3*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_x[_qp],2) + 3*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*std::pow(_polar_y[_qp],2))*_t42111212)/8. + 
   ((4*std::pow(_antiferrodis_A_x[_qp],2)*_polar_x[_qp]*_polar_y[_qp] + 4*std::pow(_antiferrodis_A_y[_qp],2)*_polar_x[_qp]*_polar_y[_qp])*_t42112211)/6. + (2*std::pow(_antiferrodis_A_z[_qp],2)*_polar_x[_qp]*_polar_y[_qp]*_t42112233)/3. + ((2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_polar_z[_qp],2))*_t42123312)/24.);
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  ((_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_t1212)/4. + ((std::pow(_antiferrodis_A_x[_qp],3)*_antiferrodis_A_z[_qp] + _antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],3))*_t24121112)/8. + (_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*_t24121233)/24. + ((3*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_x[_qp],2) + 3*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_z[_qp]*std::pow(_polar_z[_qp],2))*_t42111212)/8. + 
   ((4*std::pow(_antiferrodis_A_x[_qp],2)*_polar_x[_qp]*_polar_z[_qp] + 4*std::pow(_antiferrodis_A_z[_qp],2)*_polar_x[_qp]*_polar_z[_qp])*_t42112211)/6. + (2*std::pow(_antiferrodis_A_y[_qp],2)*_polar_x[_qp]*_polar_z[_qp]*_t42112233)/3. + ((2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_y[_qp],2) + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312)/24.);
    }
    else if (jvar == _antiferrodis_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_x[_qp]*_polar_x[_qp]*_t1111 + ((_antiferrodis_A_y[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp])*_t1212)/4. + 8*std::pow(_antiferrodis_A_x[_qp],3)*_polar_x[_qp]*_t24111111 + (2*_antiferrodis_A_x[_qp]*(std::pow(_antiferrodis_A_y[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_polar_x[_qp]*_t24111122)/3. + ((std::pow(_antiferrodis_A_y[_qp],3)*_polar_y[_qp] + std::pow(_antiferrodis_A_z[_qp],3)*_polar_z[_qp] + 3*std::pow(_antiferrodis_A_x[_qp],2)*(_antiferrodis_A_y[_qp]*_polar_y[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp]))*_t24121112)/8. + 
   ((_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_y[_qp] + std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*_polar_z[_qp])*_t24121233)/24. + 8*_antiferrodis_A_x[_qp]*std::pow(_polar_x[_qp],3)*_t42111111 + ((2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_polar_y[_qp],3) + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_z[_qp],3) + 3*std::pow(_polar_x[_qp],2)*(std::pow(_antiferrodis_A_y[_qp],2)*_polar_y[_qp] + std::pow(_antiferrodis_A_z[_qp],2)*_polar_z[_qp]))*_t42111212)/8. + 
   (2*_antiferrodis_A_x[_qp]*_polar_x[_qp]*(std::pow(_polar_y[_qp],2) + std::pow(_polar_z[_qp],2))*_t42112211)/3. + ((_antiferrodis_A_z[_qp]*std::pow(_polar_y[_qp],2)*_polar_z[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp]*std::pow(_polar_z[_qp],2))*_t42123312)/24. + 16*std::pow(_antiferrodis_A_x[_qp],3)*std::pow(_polar_x[_qp],3)*_t4411111111 + 12*_antiferrodis_A_x[_qp]*std::pow(_polar_x[_qp],5)*_t6211111111);
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_t1122 + (_antiferrodis_A_x[_qp]*_polar_y[_qp]*_t1212)/4. + (2*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_t24111122)/3. + 8*std::pow(_antiferrodis_A_y[_qp],3)*_polar_x[_qp]*_t24112222 + (2*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_x[_qp]*_t24112233)/3. + ((std::pow(_antiferrodis_A_x[_qp],3)*_polar_y[_qp] + 3*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_polar_y[_qp])*_t24121112)/8. + 
   ((_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_y[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_z[_qp])*_t24121233)/24. + 12*std::pow(_antiferrodis_A_y[_qp],5)*_polar_x[_qp]*_t2611222222 + 8*_antiferrodis_A_y[_qp]*std::pow(_polar_x[_qp],3)*_t42111122 + ((6*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_polar_x[_qp],2)*_polar_y[_qp] + std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_y[_qp],3))*_t42111212)/8. + 
   (2*_antiferrodis_A_y[_qp]*_polar_x[_qp]*std::pow(_polar_y[_qp],2)*_t42112211)/3. + (2*_antiferrodis_A_y[_qp]*_polar_x[_qp]*std::pow(_polar_z[_qp],2)*_t42112233)/3. + ((2*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiferrodis_A_x[_qp]*_polar_y[_qp]*std::pow(_polar_z[_qp],2))*_t42123312)/24. + 16*std::pow(_antiferrodis_A_y[_qp],3)*std::pow(_polar_x[_qp],3)*_t4411112222 + 12*_antiferrodis_A_y[_qp]*std::pow(_polar_x[_qp],5)*_t6211111122);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_t1122 + (_antiferrodis_A_x[_qp]*_polar_z[_qp]*_t1212)/4. + (2*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_t24111122)/3. + 8*std::pow(_antiferrodis_A_z[_qp],3)*_polar_x[_qp]*_t24112222 + (2*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_t24112233)/3. + ((std::pow(_antiferrodis_A_x[_qp],3)*_polar_z[_qp] + 3*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_z[_qp])*_t24121112)/8. + 
   ((2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp] + _antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_polar_z[_qp])*_t24121233)/24. + 12*std::pow(_antiferrodis_A_z[_qp],5)*_polar_x[_qp]*_t2611222222 + 8*_antiferrodis_A_z[_qp]*std::pow(_polar_x[_qp],3)*_t42111122 + ((6*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_x[_qp],2)*_polar_z[_qp] + std::pow(_antiferrodis_A_x[_qp],2)*std::pow(_polar_z[_qp],3))*_t42111212)/8. + 
   (2*_antiferrodis_A_z[_qp]*_polar_x[_qp]*std::pow(_polar_z[_qp],2)*_t42112211)/3. + (2*_antiferrodis_A_z[_qp]*_polar_x[_qp]*std::pow(_polar_y[_qp],2)*_t42112233)/3. + ((2*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiferrodis_A_x[_qp]*std::pow(_polar_y[_qp],2)*_polar_z[_qp])*_t42123312)/24. + 16*std::pow(_antiferrodis_A_z[_qp],3)*std::pow(_polar_x[_qp],3)*_t4411112222 + 12*_antiferrodis_A_z[_qp]*std::pow(_polar_x[_qp],5)*_t6211111122);
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
      return _test[_i][_qp] * _phi[_j][_qp] *  ((_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_t1212)/4. + ((std::pow(_antiferrodis_A_x[_qp],3)*_antiferrodis_A_y[_qp] + _antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],3))*_t24121112)/8. + (_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_t24121233)/24. + ((3*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_x[_qp],2) + 3*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*std::pow(_polar_y[_qp],2))*_t42111212)/8. + 
   ((4*std::pow(_antiferrodis_A_x[_qp],2)*_polar_x[_qp]*_polar_y[_qp] + 4*std::pow(_antiferrodis_A_y[_qp],2)*_polar_x[_qp]*_polar_y[_qp])*_t42112211)/6. + (2*std::pow(_antiferrodis_A_z[_qp],2)*_polar_x[_qp]*_polar_y[_qp]*_t42112233)/3. + ((2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_z[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_polar_z[_qp],2))*_t42123312)/24.);
    }
    else if (jvar == _polar_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  ((_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_t1212)/4. + ((std::pow(_antiferrodis_A_y[_qp],3)*_antiferrodis_A_z[_qp] + _antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],3))*_t24121112)/8. + (std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_t24121233)/24. + ((3*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_y[_qp],2) + 3*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*std::pow(_polar_z[_qp],2))*_t42111212)/8. + 
   ((4*std::pow(_antiferrodis_A_y[_qp],2)*_polar_y[_qp]*_polar_z[_qp] + 4*std::pow(_antiferrodis_A_z[_qp],2)*_polar_y[_qp]*_polar_z[_qp])*_t42112211)/6. + (2*std::pow(_antiferrodis_A_x[_qp],2)*_polar_y[_qp]*_polar_z[_qp]*_t42112233)/3. + ((_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_x[_qp],2) + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_z[_qp])*_t42123312)/24.);
    }
    else if (jvar == _antiferrodis_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_x[_qp]*_polar_y[_qp]*_t1122 + (_antiferrodis_A_y[_qp]*_polar_x[_qp]*_t1212)/4. + (2*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_polar_y[_qp]*_t24111122)/3. + 8*std::pow(_antiferrodis_A_x[_qp],3)*_polar_y[_qp]*_t24112222 + (2*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_y[_qp]*_t24112233)/3. + ((3*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_polar_x[_qp] + std::pow(_antiferrodis_A_y[_qp],3)*_polar_x[_qp])*_t24121112)/8. + 
   ((_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_x[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_z[_qp])*_t24121233)/24. + 12*std::pow(_antiferrodis_A_x[_qp],5)*_polar_y[_qp]*_t2611222222 + 8*_antiferrodis_A_x[_qp]*std::pow(_polar_y[_qp],3)*_t42111122 + ((std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_x[_qp],3) + 6*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_x[_qp]*std::pow(_polar_y[_qp],2))*_t42111212)/8. + 
   (2*_antiferrodis_A_x[_qp]*std::pow(_polar_x[_qp],2)*_polar_y[_qp]*_t42112211)/3. + (2*_antiferrodis_A_x[_qp]*_polar_y[_qp]*std::pow(_polar_z[_qp],2)*_t42112233)/3. + ((2*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp] + _antiferrodis_A_y[_qp]*_polar_x[_qp]*std::pow(_polar_z[_qp],2))*_t42123312)/24. + 8*std::pow(_antiferrodis_A_x[_qp],3)*_polar_y[_qp]*_t4411112222 + 12*_antiferrodis_A_x[_qp]*std::pow(_polar_y[_qp],5)*_t6211111122);
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_y[_qp]*_polar_y[_qp]*_t1111 + ((_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp])*_t1212)/4. + 8*std::pow(_antiferrodis_A_y[_qp],3)*_polar_y[_qp]*_t24111111 + (2*_antiferrodis_A_y[_qp]*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_z[_qp],2))*_polar_y[_qp]*_t24111122)/3. + ((std::pow(_antiferrodis_A_x[_qp],3)*_polar_x[_qp] + std::pow(_antiferrodis_A_z[_qp],3)*_polar_z[_qp] + 3*std::pow(_antiferrodis_A_y[_qp],2)*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_z[_qp]*_polar_z[_qp]))*_t24121112)/8. + 
   ((_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_x[_qp] + std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_z[_qp]*_polar_z[_qp])*_t24121233)/24. + 8*_antiferrodis_A_y[_qp]*std::pow(_polar_y[_qp],3)*_t42111111 + ((2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*std::pow(_polar_x[_qp],3) + 2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_z[_qp],3) + 3*std::pow(_polar_y[_qp],2)*(std::pow(_antiferrodis_A_x[_qp],2)*_polar_x[_qp] + std::pow(_antiferrodis_A_z[_qp],2)*_polar_z[_qp]))*_t42111212)/8. + 
   (2*_antiferrodis_A_y[_qp]*_polar_y[_qp]*(std::pow(_polar_x[_qp],2) + std::pow(_polar_z[_qp],2))*_t42112211)/3. + ((_antiferrodis_A_z[_qp]*std::pow(_polar_x[_qp],2)*_polar_z[_qp] + _antiferrodis_A_x[_qp]*_polar_x[_qp]*std::pow(_polar_z[_qp],2))*_t42123312)/24. + 16*std::pow(_antiferrodis_A_y[_qp],3)*std::pow(_polar_y[_qp],3)*_t4411111111 + 12*_antiferrodis_A_y[_qp]*std::pow(_polar_y[_qp],5)*_t6211111111);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_t1122 + (_antiferrodis_A_y[_qp]*_polar_z[_qp]*_t1212)/4. + (2*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_t24111122)/3. + 8*std::pow(_antiferrodis_A_z[_qp],3)*_polar_y[_qp]*_t24112222 + (2*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_z[_qp]*_polar_y[_qp]*_t24112233)/3. + ((std::pow(_antiferrodis_A_y[_qp],3)*_polar_z[_qp] + 3*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_z[_qp])*_t24121112)/8. + 
   ((2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp] + std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_polar_z[_qp])*_t24121233)/24. + 12*std::pow(_antiferrodis_A_z[_qp],5)*_polar_y[_qp]*_t2611222222 + 8*_antiferrodis_A_z[_qp]*std::pow(_polar_y[_qp],3)*_t42111122 + ((6*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_y[_qp],2)*_polar_z[_qp] + std::pow(_antiferrodis_A_y[_qp],2)*std::pow(_polar_z[_qp],3))*_t42111212)/8. + 
   (2*_antiferrodis_A_z[_qp]*_polar_y[_qp]*std::pow(_polar_z[_qp],2)*_t42112211)/3. + (2*_antiferrodis_A_z[_qp]*std::pow(_polar_x[_qp],2)*_polar_y[_qp]*_t42112233)/3. + ((_antiferrodis_A_y[_qp]*std::pow(_polar_x[_qp],2)*_polar_z[_qp] + 2*_antiferrodis_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312)/24. + 8*std::pow(_antiferrodis_A_z[_qp],3)*_polar_y[_qp]*_t4411112222 + 12*_antiferrodis_A_z[_qp]*std::pow(_polar_y[_qp],5)*_t6211111122);
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
      return _test[_i][_qp] * _phi[_j][_qp] *  ((_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_t1212)/4. + ((std::pow(_antiferrodis_A_x[_qp],3)*_antiferrodis_A_z[_qp] + _antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],3))*_t24121112)/8. + (_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*_t24121233)/24. + ((3*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_x[_qp],2) + 3*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_z[_qp]*std::pow(_polar_z[_qp],2))*_t42111212)/8. + 
   ((4*std::pow(_antiferrodis_A_x[_qp],2)*_polar_x[_qp]*_polar_z[_qp] + 4*std::pow(_antiferrodis_A_z[_qp],2)*_polar_x[_qp]*_polar_z[_qp])*_t42112211)/6. + (2*std::pow(_antiferrodis_A_y[_qp],2)*_polar_x[_qp]*_polar_z[_qp]*_t42112233)/3. + ((2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + _antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_y[_qp],2) + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312)/24.);
    }
    else if (jvar == _polar_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  ((_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_t1212)/4. + ((std::pow(_antiferrodis_A_y[_qp],3)*_antiferrodis_A_z[_qp] + _antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],3))*_t24121112)/8. + (std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_t24121233)/24. + ((3*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_y[_qp],2) + 3*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*std::pow(_polar_z[_qp],2))*_t42111212)/8. + 
   ((4*std::pow(_antiferrodis_A_y[_qp],2)*_polar_y[_qp]*_polar_z[_qp] + 4*std::pow(_antiferrodis_A_z[_qp],2)*_polar_y[_qp]*_polar_z[_qp])*_t42112211)/6. + (2*std::pow(_antiferrodis_A_x[_qp],2)*_polar_y[_qp]*_polar_z[_qp]*_t42112233)/3. + ((_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_x[_qp],2) + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*_polar_y[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_z[_qp])*_t42123312)/24.);
    }
    else if (jvar == _antiferrodis_A_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_x[_qp]*_polar_z[_qp]*_t1122 + (_antiferrodis_A_z[_qp]*_polar_x[_qp]*_t1212)/4. + (2*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_z[_qp]*_t24111122)/3. + 8*std::pow(_antiferrodis_A_x[_qp],3)*_polar_z[_qp]*_t24112222 + (2*_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_polar_z[_qp]*_t24112233)/3. + ((3*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_z[_qp]*_polar_x[_qp] + std::pow(_antiferrodis_A_z[_qp],3)*_polar_x[_qp])*_t24121112)/8. + ((std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*_polar_x[_qp] + 2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp])*_t24121233)/24. + 12*std::pow(_antiferrodis_A_x[_qp],5)*_polar_z[_qp]*_t2611222222 + 8*_antiferrodis_A_x[_qp]*std::pow(_polar_z[_qp],3)*_t42111122 + ((std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_x[_qp],3) + 6*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp]*std::pow(_polar_z[_qp],2))*_t42111212)/8. + (2*_antiferrodis_A_x[_qp]*std::pow(_polar_x[_qp],2)*_polar_z[_qp]*_t42112211)/3. + (2*_antiferrodis_A_x[_qp]*std::pow(_polar_y[_qp],2)*_polar_z[_qp]*_t42112233)/3. + ((_antiferrodis_A_z[_qp]*_polar_x[_qp]*std::pow(_polar_y[_qp],2) + 2*_antiferrodis_A_y[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312)/24. + 16*std::pow(_antiferrodis_A_x[_qp],3)*std::pow(_polar_z[_qp],3)*_t4411112222 + 12*_antiferrodis_A_x[_qp]*std::pow(_polar_z[_qp],5)*_t6211111122);
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_y[_qp]*_polar_z[_qp]*_t1122 + (_antiferrodis_A_z[_qp]*_polar_y[_qp]*_t1212)/4. + (2*_antiferrodis_A_y[_qp]*std::pow(_antiferrodis_A_z[_qp],2)*_polar_z[_qp]*_t24111122)/3. + 8*std::pow(_antiferrodis_A_y[_qp],3)*_polar_z[_qp]*_t24112222 + (2*std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_polar_z[_qp]*_t24112233)/3. + ((3*std::pow(_antiferrodis_A_y[_qp],2)*_antiferrodis_A_z[_qp]*_polar_y[_qp] + std::pow(_antiferrodis_A_z[_qp],3)*_polar_y[_qp])*_t24121112)/8. + ((2*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_x[_qp] + std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_z[_qp]*_polar_y[_qp])*_t24121233)/24. + 12*std::pow(_antiferrodis_A_y[_qp],5)*_polar_z[_qp]*_t2611222222 + 8*_antiferrodis_A_y[_qp]*std::pow(_polar_z[_qp],3)*_t42111122 + ((std::pow(_antiferrodis_A_z[_qp],2)*std::pow(_polar_y[_qp],3) + 6*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_polar_y[_qp]*std::pow(_polar_z[_qp],2))*_t42111212)/8. + (2*_antiferrodis_A_y[_qp]*std::pow(_polar_y[_qp],2)*_polar_z[_qp]*_t42112211)/3. + (2*_antiferrodis_A_y[_qp]*std::pow(_polar_x[_qp],2)*_polar_z[_qp]*_t42112233)/3. + ((_antiferrodis_A_z[_qp]*std::pow(_polar_x[_qp],2)*_polar_y[_qp] + 2*_antiferrodis_A_x[_qp]*_polar_x[_qp]*_polar_y[_qp]*_polar_z[_qp])*_t42123312)/24. + 16*std::pow(_antiferrodis_A_y[_qp],3)*std::pow(_polar_z[_qp],3)*_t4411112222 + 12*_antiferrodis_A_y[_qp]*std::pow(_polar_z[_qp],5)*_t6211111122);
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (4*_antiferrodis_A_z[_qp]*_polar_z[_qp]*_t1111 + ((_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp])*_t1212)/4. + 8*std::pow(_antiferrodis_A_z[_qp],3)*_polar_z[_qp]*_t24111111 + (2*(std::pow(_antiferrodis_A_x[_qp],2) + std::pow(_antiferrodis_A_y[_qp],2))*_antiferrodis_A_z[_qp]*_polar_z[_qp]*_t24111122)/3. + ((std::pow(_antiferrodis_A_x[_qp],3)*_polar_x[_qp] + std::pow(_antiferrodis_A_y[_qp],3)*_polar_y[_qp] + 3*std::pow(_antiferrodis_A_z[_qp],2)*(_antiferrodis_A_x[_qp]*_polar_x[_qp] + _antiferrodis_A_y[_qp]*_polar_y[_qp]))*_t24121112)/8. + ((_antiferrodis_A_x[_qp]*std::pow(_antiferrodis_A_y[_qp],2)*_polar_x[_qp] + std::pow(_antiferrodis_A_x[_qp],2)*_antiferrodis_A_y[_qp]*_polar_y[_qp])*_t24121233)/24. + 8*_antiferrodis_A_z[_qp]*std::pow(_polar_z[_qp],3)*_t42111111 + ((2*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_x[_qp],3) + 2*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*std::pow(_polar_y[_qp],3) + 3*(std::pow(_antiferrodis_A_x[_qp],2)*_polar_x[_qp] + std::pow(_antiferrodis_A_y[_qp],2)*_polar_y[_qp])*std::pow(_polar_z[_qp],2))*_t42111212)/8. + (2*_antiferrodis_A_z[_qp]*(std::pow(_polar_x[_qp],2) + std::pow(_polar_y[_qp],2))*_polar_z[_qp]*_t42112211)/3. + ((_antiferrodis_A_y[_qp]*std::pow(_polar_x[_qp],2)*_polar_y[_qp] + _antiferrodis_A_x[_qp]*_polar_x[_qp]*std::pow(_polar_y[_qp],2))*_t42123312)/24. + 16*std::pow(_antiferrodis_A_z[_qp],3)*std::pow(_polar_z[_qp],3)*_t4411111111 + 12*_antiferrodis_A_z[_qp]*std::pow(_polar_z[_qp],5)*_t6211111111);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
