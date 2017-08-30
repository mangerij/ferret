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

#include "PiezoelectricStrainCharge.h"
#include "ComputePiezoTensor.h"
#include "PiezostrictiveTensorTools.h"

class PiezoelectricStrainCharge;

template<>
InputParameters validParams<PiezoelectricStrainCharge>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("disp_x", "The x component of the displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the displacement");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

PiezoelectricStrainCharge::PiezoelectricStrainCharge(const InputParameters & parameters)
  :Kernel(parameters),
   _piezo_tensor(getMaterialProperty<RankThreeTensor>("piezo_tensor")),
   _piezostrictive_tensor(getMaterialProperty<RankThreeTensor>("piezostrictive_tensor")),
   _disp_x_var(coupled("disp_x")),
   _disp_y_var(coupled("disp_y")),
   _disp_z_var(coupled("disp_z")),
   _disp_x_grad(coupledGradient("disp_x")),
   _disp_y_grad(coupledGradient("disp_y")),
   _disp_z_grad(coupledGradient("disp_z")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
PiezoelectricStrainCharge::computeQpResidual()
{
  Real sum = 0.0;
  for(unsigned int j = 0; j < 3; ++j)
  {
    sum += - _grad_test[_i][_qp](0) * (_piezostrictive_tensor[_qp](0,0,j) * _disp_x_grad[_qp](j) + _piezostrictive_tensor[_qp](0,1,j) * _disp_y_grad[_qp](j) + _piezostrictive_tensor[_qp](0,2,j) * _disp_z_grad[_qp](j)) + _grad_test[_i][_qp](1) * (_piezostrictive_tensor[_qp](1,0,j) * _disp_x_grad[_qp](j) + _piezostrictive_tensor[_qp](1,1,j) * _disp_y_grad[_qp](j) + _piezostrictive_tensor[_qp](1,2,j) * _disp_z_grad[_qp](j)) + _grad_test[_i][_qp](2) * (_piezostrictive_tensor[_qp](2,0,j) * _disp_x_grad[_qp](j) + _piezostrictive_tensor[_qp](2,1,j) * _disp_y_grad[_qp](j) + _piezostrictive_tensor[_qp](2,2,j) * _disp_z_grad[_qp](j));
  }
  return sum;
}


Real
PiezoelectricStrainCharge::computeQpJacobian()
{
  return 0.0;
}

Real
PiezoelectricStrainCharge::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real sum = 0.0;
  if(jvar == _disp_x_var || jvar == _disp_y_var || jvar == _disp_z_var)
  {
    if (jvar == _disp_x_var)
      {
        for(unsigned int j = 0; j < 3; ++j)
          {
            sum += - _grad_test[_i][_qp](0) * (_piezostrictive_tensor[_qp](0,0,j) * _grad_phi[_j][_qp](j) ) + _grad_test[_i][_qp](1) * (_piezostrictive_tensor[_qp](1,0,j) * _grad_phi[_j][_qp](j)) + _grad_test[_i][_qp](2) * (_piezostrictive_tensor[_qp](2,0,j) * _grad_phi[_j][_qp](j));
          }
      }
    else if (jvar == _disp_y_var)
      {
        for(unsigned int j = 0; j < 3; ++j)
          {
            sum += - _grad_test[_i][_qp](0) * ( _piezostrictive_tensor[_qp](0,1,j) * _grad_phi[_j][_qp](j) ) + _grad_test[_i][_qp](1) * ( _piezostrictive_tensor[_qp](1,1,j) * _grad_phi[_j][_qp](j) ) + _grad_test[_i][_qp](2) * (_piezostrictive_tensor[_qp](2,1,j) * _grad_phi[_j][_qp](j));
          }
      }
    else if (jvar == _disp_z_var)
      {
        for(unsigned int j = 0; j < 3; ++j)
          {
            sum += - _grad_test[_i][_qp](0) * (_piezostrictive_tensor[_qp](0,2,j) * _grad_phi[_j][_qp](j)) + _grad_test[_i][_qp](1) * (_piezostrictive_tensor[_qp](1,2,j) * _grad_phi[_j][_qp](j)) + _grad_test[_i][_qp](2) * (_piezostrictive_tensor[_qp](2,2,j) * _grad_phi[_j][_qp](j));
          }
      }
    return sum;
  }
  else
  {
    return 0.0;
  }
}

