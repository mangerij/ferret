/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#include "RotatedBulkEnergyDerivative.h"
#include<cmath>

template<>
InputParameters validParams<RotatedBulkEnergyDerivative>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("euler_angle_3", 0.0, "Euler angle in direction 3");
  params.addRequiredParam<Real>("alpha1", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha11", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha12", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha123", "The coefficients of the Landau expansion");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

RotatedBulkEnergyDerivative::RotatedBulkEnergyDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _Euler_angles(getParam<Real>("euler_angle_1"),
                 getParam<Real>("euler_angle_2"),
                 getParam<Real>("euler_angle_3")),
   _alpha1(getParam<Real>("alpha1")),
   _alpha11(getParam<Real>("alpha11")),
   _alpha12(getParam<Real>("alpha12")),
   _alpha111(getParam<Real>("alpha111")),
   _alpha112(getParam<Real>("alpha112")),
   _alpha123(getParam<Real>("alpha123")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
RotatedBulkEnergyDerivative::computeQpResidual()
{
  
  if(_component == 0)
    {
      RotationTensor R(_Euler_angles);
      return _test[_i][_qp] * ((2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)))*_alpha1 +
   (4*R(0,0)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3) + 4*R(1,0)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3) + 4*R(2,0)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3))*_alpha11 +
   (6*R(0,0)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),5) + 6*R(1,0)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),5) + 6*R(2,0)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),5))*_alpha111 +
   (4*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2))*R(2,0)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) +
      (2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4) +
      std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4)*(2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) +
      std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4)*(2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) +
      4*R(1,0)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)) +
      4*R(0,0)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)))*_alpha112 +
   (2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) +
      2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 2*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) +
      2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha12 +
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 +
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 +
   2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123);
    }
  else if(_component == 1)
    {
      RotationTensor R(_Euler_angles);
      return _test[_i][_qp] * ((2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)))*_alpha1 +
   (4*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3) + 4*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3) + 4*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3))*_alpha11 +
   (6*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),5) + 6*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),5) + 6*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),5))*_alpha111 +
   (4*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2))*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) +
      (2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4) +
      std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4)*(2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) +
      std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4)*(2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) +
      4*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)) +
      4*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)))*_alpha112 +
   (2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) +
      2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 2*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) +
      2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha12 +
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 +
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 +
   2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123);
    }
  else
    {
      RotationTensor R(_Euler_angles);
      return _test[_i][_qp] * ((2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)))*_alpha1 +
   (4*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3) + 4*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3) + 4*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3))*_alpha11 +
   (6*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),5) + 6*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),5) + 6*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),5))*_alpha111 +
   (4*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2))*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) +
      (2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4) +
      std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4)*(2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) +
      std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4)*(2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) +
      4*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)) +
      4*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)))*_alpha112 +
   (2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) +
      2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 2*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) +
      2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha12 +
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 +
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 +
   2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123);
    }
}

Real
RotatedBulkEnergyDerivative::computeQpJacobian()
{
  if(_component == 0)
    {
      RotationTensor R(_Euler_angles);
      return _test[_i][_qp] * _phi[_j][_qp] * ((2*std::pow(R(0,0),2) + 2*std::pow(R(1,0),2) + 2*std::pow(R(2,0),2))*_alpha1 + (12*std::pow(R(0,0),2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + 12*std::pow(R(1,0),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) +
      12*std::pow(R(2,0),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha11 + (30*std::pow(R(0,0),2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4) + 30*std::pow(R(1,0),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4) +
      30*std::pow(R(2,0),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4))*_alpha111 + (std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4)*(2*std::pow(R(0,0),2) + 2*std::pow(R(2,0),2)) +
      std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4)*(2*std::pow(R(1,0),2) + 2*std::pow(R(2,0),2)) + 12*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2))*std::pow(R(2,0),2)*
       std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 8*(2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,0)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) +
      (2*std::pow(R(0,0),2) + 2*std::pow(R(1,0),2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4) + 8*R(1,0)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) +
      8*R(0,0)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) +
      12*std::pow(R(1,0),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)) +
      12*std::pow(R(0,0),2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*(std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)))*_alpha112 +
   (2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(R(1,0),2) + 8*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*std::pow(R(0,0),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) +
      2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(R(2,0),2) + 2*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(R(2,0),2) + 8*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) +
      8*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 2*std::pow(R(0,0),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 2*std::pow(R(1,0),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha12 +
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(R(2,0),2)*_alpha123 +
   8*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 +
   8*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 +
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(R(1,0),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 +
   8*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 +
   2*std::pow(R(0,0),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123);
    }
  else if(_component == 1)
    {
      RotationTensor R(_Euler_angles);
      return _test[_i][_qp] * _phi[_j][_qp] * ((2*std::pow(R(0,1),2) + 2*std::pow(R(1,1),2) + 2*std::pow(R(2,1),2))*_alpha1 + (12*std::pow(R(0,1),2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + 12*std::pow(R(1,1),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) +
      12*std::pow(R(2,1),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha11 + (30*std::pow(R(0,1),2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4) + 30*std::pow(R(1,1),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4) +
      30*std::pow(R(2,1),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4))*_alpha111 + (std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4)*(2*std::pow(R(0,1),2) + 2*std::pow(R(2,1),2)) +
      std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4)*(2*std::pow(R(1,1),2) + 2*std::pow(R(2,1),2)) + 12*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2))*std::pow(R(2,1),2)*
       std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 8*(2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) +
      (2*std::pow(R(0,1),2) + 2*std::pow(R(1,1),2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4) + 8*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) +
      8*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) +
      12*std::pow(R(1,1),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)) +
      12*std::pow(R(0,1),2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*(std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)))*_alpha112 +
   (2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(R(1,1),2) + 8*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*std::pow(R(0,1),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) +
      2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(R(2,1),2) + 2*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(R(2,1),2) + 8*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) +
      8*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 2*std::pow(R(0,1),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 2*std::pow(R(1,1),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha12 +
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(R(2,1),2)*_alpha123 +
   8*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 +
   8*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 +
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(R(1,1),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 +
   8*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 +
   2*std::pow(R(0,1),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123);
    }
  else
    {
      RotationTensor R(_Euler_angles);
      return _test[_i][_qp] * _phi[_j][_qp] * ((2*std::pow(R(0,2),2) + 2*std::pow(R(1,2),2) + 2*std::pow(R(2,2),2))*_alpha1 + (12*std::pow(R(0,2),2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + 12*std::pow(R(1,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) +
      12*std::pow(R(2,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha11 + (30*std::pow(R(0,2),2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4) + 30*std::pow(R(1,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4) +
      30*std::pow(R(2,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4))*_alpha111 + (12*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2))*std::pow(R(2,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) +
      8*(2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + (2*std::pow(R(0,2),2) + 2*std::pow(R(1,2),2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4) +
      std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4)*(2*std::pow(R(0,2),2) + 2*std::pow(R(2,2),2)) + std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4)*(2*std::pow(R(1,2),2) + 2*std::pow(R(2,2),2)) +
      8*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) +
      8*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) +
      12*std::pow(R(1,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)) +
      12*std::pow(R(0,2),2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*(std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)))*_alpha112 +
   (2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(R(1,2),2) + 8*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*std::pow(R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) +
      2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(R(2,2),2) + 2*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(R(2,2),2) + 8*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) +
      8*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 2*std::pow(R(0,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 2*std::pow(R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha12 +
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(R(2,2),2)*_alpha123 +
   8*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 +
   8*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 +
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 +
   8*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 +
   2*std::pow(R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123);
    }
}

Real
RotatedBulkEnergyDerivative::computeQpOffDiagJacobian(unsigned int jvar)

{
  if (jvar == _polar_x_var)
  {
    if (_component == 1)
    {
      RotationTensor R(_Euler_angles);
      return _test[_i][_qp] * _phi[_j][_qp] * ((2*R(0,0)*R(0,1) + 2*R(1,0)*R(1,1) + 2*R(2,0)*R(2,1))*_alpha1 + (12*R(0,0)*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + 12*R(1,0)*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + 12*R(2,0)*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha11 + 
   (30*R(0,0)*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4) + 30*R(1,0)*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4) + 30*R(2,0)*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4))*_alpha111 + 
   (std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4)*(2*R(0,0)*R(0,1) + 2*R(2,0)*R(2,1)) + std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4)*(2*R(1,0)*R(1,1) + 2*R(2,0)*R(2,1)) + 
      12*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2))*R(2,0)*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 
      4*(2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,0)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + 
      4*(2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + (2*R(0,0)*R(0,1) + 2*R(1,0)*R(1,1))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4) + 
      4*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(1,0)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(0,0)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      12*R(1,0)*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)) + 
      12*R(0,0)*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*(std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)))*_alpha112 + 
   (2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*R(1,1) + 4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 
      2*R(0,0)*R(0,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + 2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(2,0)*R(2,1) + 2*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*R(2,1) + 
      4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 4*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 
      4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 4*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 2*R(0,0)*R(0,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 
      2*R(1,0)*R(1,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha12 + 2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*R(2,1)*_alpha123 + 
   4*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*R(1,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 
   4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 
   4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 2*R(0,0)*R(0,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123);
    }
    else if (_component == 2)
    {
      RotationTensor R(_Euler_angles);
      return _test[_i][_qp] * _phi[_j][_qp] * ((2*R(0,0)*R(0,2) + 2*R(1,0)*R(1,2) + 2*R(2,0)*R(2,2))*_alpha1 + (12*R(0,0)*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + 12*R(1,0)*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + 12*R(2,0)*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha11 + 
   (30*R(0,0)*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4) + 30*R(1,0)*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4) + 30*R(2,0)*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4))*_alpha111 + 
   (12*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2))*R(2,0)*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 
      4*(2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,0)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + 
      4*(2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + (2*R(0,0)*R(0,2) + 2*R(1,0)*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4) + 
      std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4)*(2*R(0,0)*R(0,2) + 2*R(2,0)*R(2,2)) + std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4)*(2*R(1,0)*R(1,2) + 2*R(2,0)*R(2,2)) + 
      4*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(1,0)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(0,0)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      12*R(1,0)*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)) + 
      12*R(0,0)*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*(std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)))*_alpha112 + 
   (2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*R(1,2) + 4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 
      2*R(0,0)*R(0,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + 2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(2,0)*R(2,2) + 2*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*R(2,2) + 
      4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 4*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 
      4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 4*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 2*R(0,0)*R(0,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 
      2*R(1,0)*R(1,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha12 + 2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*R(2,2)*_alpha123 + 
   4*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*R(1,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 
   4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 
   4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 2*R(0,0)*R(0,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123);
    }
    else
    {
      return 0.0;
    }
  }
  else if (jvar == _polar_y_var)
  {
    if (_component == 0)
    {
      RotationTensor R(_Euler_angles);
      return _test[_i][_qp] * _phi[_j][_qp] * ((2*R(0,0)*R(0,1) + 2*R(1,0)*R(1,1) + 2*R(2,0)*R(2,1))*_alpha1 + (12*R(0,0)*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + 12*R(1,0)*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + 12*R(2,0)*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha11 + 
   (30*R(0,0)*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4) + 30*R(1,0)*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4) + 30*R(2,0)*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4))*_alpha111 + 
   (std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4)*(2*R(0,0)*R(0,1) + 2*R(2,0)*R(2,1)) + std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4)*(2*R(1,0)*R(1,1) + 2*R(2,0)*R(2,1)) + 
      12*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2))*R(2,0)*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 
      4*(2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,0)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + 
      4*(2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + (2*R(0,0)*R(0,1) + 2*R(1,0)*R(1,1))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4) + 
      4*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(1,0)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(0,0)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      12*R(1,0)*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)) + 
      12*R(0,0)*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*(std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)))*_alpha112 + 
   (2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*R(1,1) + 4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 
      2*R(0,0)*R(0,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + 2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(2,0)*R(2,1) + 2*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*R(2,1) + 
      4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 4*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 
      4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 4*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 2*R(0,0)*R(0,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 
      2*R(1,0)*R(1,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha12 + 2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*R(2,1)*_alpha123 + 
   4*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*R(1,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 
   4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 
   4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 2*R(0,0)*R(0,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123);
    }
    else if (_component == 2)
    {
      RotationTensor R(_Euler_angles);
      return _test[_i][_qp] * _phi[_j][_qp] * ((2*R(0,1)*R(0,2) + 2*R(1,1)*R(1,2) + 2*R(2,1)*R(2,2))*_alpha1 + (12*R(0,1)*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + 12*R(1,1)*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + 12*R(2,1)*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha11 + 
   (30*R(0,1)*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4) + 30*R(1,1)*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4) + 30*R(2,1)*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4))*_alpha111 + 
   (12*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2))*R(2,1)*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 
      4*(2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + 
      4*(2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + (2*R(0,1)*R(0,2) + 2*R(1,1)*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4) + 
      std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4)*(2*R(0,1)*R(0,2) + 2*R(2,1)*R(2,2)) + std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4)*(2*R(1,1)*R(1,2) + 2*R(2,1)*R(2,2)) + 
      4*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      12*R(1,1)*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)) + 
      12*R(0,1)*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*(std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)))*_alpha112 + 
   (2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,1)*R(1,2) + 4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 
      2*R(0,1)*R(0,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + 2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(2,1)*R(2,2) + 2*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,1)*R(2,2) + 
      4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 4*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 
      4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 4*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 2*R(0,1)*R(0,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 
      2*R(1,1)*R(1,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha12 + 2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,1)*R(2,2)*_alpha123 + 
   4*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,1)*R(1,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 
   4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 
   4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 2*R(0,1)*R(0,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123);
    }
    else
    {
      return 0.0;
    }
  }
  else if (jvar == _polar_z_var)
  {
    if (_component == 0)
    {
      RotationTensor R(_Euler_angles);
      return _test[_i][_qp] * _phi[_j][_qp] * ((2*R(0,0)*R(0,2) + 2*R(1,0)*R(1,2) + 2*R(2,0)*R(2,2))*_alpha1 + (12*R(0,0)*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + 12*R(1,0)*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + 12*R(2,0)*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha11 + 
   (30*R(0,0)*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4) + 30*R(1,0)*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4) + 30*R(2,0)*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4))*_alpha111 + 
   (12*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2))*R(2,0)*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 
      4*(2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,0)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + 
      4*(2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + (2*R(0,0)*R(0,2) + 2*R(1,0)*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4) + 
      std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4)*(2*R(0,0)*R(0,2) + 2*R(2,0)*R(2,2)) + std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4)*(2*R(1,0)*R(1,2) + 2*R(2,0)*R(2,2)) + 
      4*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(1,0)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(0,0)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      12*R(1,0)*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)) + 
      12*R(0,0)*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*(std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)))*_alpha112 + 
   (2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*R(1,2) + 4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 
      2*R(0,0)*R(0,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + 2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(2,0)*R(2,2) + 2*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*R(2,2) + 
      4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 4*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 
      4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 4*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 2*R(0,0)*R(0,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 
      2*R(1,0)*R(1,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha12 + 2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*R(2,2)*_alpha123 + 
   4*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,0)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,0)*R(1,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 
   4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,0)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 
   4*R(0,0)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 2*R(0,0)*R(0,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123);
    }
    else if (_component == 1)
    {
      RotationTensor R(_Euler_angles);
      return _test[_i][_qp] * _phi[_j][_qp] * ((2*R(0,1)*R(0,2) + 2*R(1,1)*R(1,2) + 2*R(2,1)*R(2,2))*_alpha1 + (12*R(0,1)*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + 12*R(1,1)*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + 12*R(2,1)*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha11 + 
   (30*R(0,1)*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4) + 30*R(1,1)*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4) + 30*R(2,1)*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4))*_alpha111 + 
   (12*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2))*R(2,1)*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 
      4*(2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,1)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + 
      4*(2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)))*R(2,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),3) + (2*R(0,1)*R(0,2) + 2*R(1,1)*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),4) + 
      std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),4)*(2*R(0,1)*R(0,2) + 2*R(2,1)*R(2,2)) + std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),4)*(2*R(1,1)*R(1,2) + 2*R(2,1)*R(2,2)) + 
      4*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(1,1)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),3)*(2*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      4*R(0,1)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),3)*(2*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 2*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))) + 
      12*R(1,1)*R(1,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*(std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)) + 
      12*R(0,1)*R(0,2)*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*(std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)))*_alpha112 + 
   (2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,1)*R(1,2) + 4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2)) + 
      2*R(0,1)*R(0,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2) + 2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(2,1)*R(2,2) + 2*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,1)*R(2,2) + 
      4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 4*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 
      4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 4*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2)) + 2*R(0,1)*R(0,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2) + 
      2*R(1,1)*R(1,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2))*_alpha12 + 2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,1)*R(2,2)*_alpha123 + 
   4*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,1)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*R(2,2)*(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2))*_alpha123 + 
   2*std::pow(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2),2)*R(1,1)*R(1,2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 
   4*R(0,2)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,1)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 
   4*R(0,1)*(_polar_x[_qp]*R(0,0) + _polar_y[_qp]*R(0,1) + _polar_z[_qp]*R(0,2))*R(1,2)*(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2))*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123 + 2*R(0,1)*R(0,2)*std::pow(_polar_x[_qp]*R(1,0) + _polar_y[_qp]*R(1,1) + _polar_z[_qp]*R(1,2),2)*std::pow(_polar_x[_qp]*R(2,0) + _polar_y[_qp]*R(2,1) + _polar_z[_qp]*R(2,2),2)*_alpha123);
    }
    else
    {
      return 0.0;
    }
  }
  else
  {
    return 0.0;
  }
}
