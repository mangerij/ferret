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

#include "ConversePiezoelectricStrain.h"
#include "ComputePiezoTensor.h"
#include "PiezostrictiveTensorTools.h"

class ConversePiezoelectricStrain;

template<>
InputParameters validParams<ConversePiezoelectricStrain>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("disp_x", "The x component of the displacement");
  params.addRequiredCoupledVar("disp_y", "The y component of the displacement");
  params.addCoupledVar("disp_z", 0.0, "The z component of the displacement");
  params.addRequiredCoupledVar("potential_int", "The electrostatic potential");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

ConversePiezoelectricStrain::ConversePiezoelectricStrain(const InputParameters & parameters)
  :Kernel(parameters),
   _piezo_tensor(getMaterialProperty<RankThreeTensor>("piezo_tensor")),
   _piezostrictive_tensor(getMaterialProperty<RankThreeTensor>("piezostrictive_tensor")),
   _component(getParam<unsigned int>("component")),
   _disp_x_var(coupled("disp_x")),
   _disp_y_var(coupled("disp_y")),
   _disp_z_var(coupled("disp_z")),
   _potential_int_var(coupled("potential_int")),
   _potential_int(coupledValue("potential_int")),
   _potential_int_second(coupledSecond("potential_int")),
   _second_phi(secondPhi()),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
ConversePiezoelectricStrain::computeQpResidual()
{
  return _test[_i][_qp] * (_piezostrictive_tensor[_qp](_component,0,0) * _potential_int_second[_qp](0,0) + _piezostrictive_tensor[_qp](_component,0,1) * _potential_int_second[_qp](0,1) + _piezostrictive_tensor[_qp](_component,0,2) * _potential_int_second[_qp](0,2) + _piezostrictive_tensor[_qp](_component,1,0) * _potential_int_second[_qp](1,0) + _piezostrictive_tensor[_qp](_component,1,1) * _potential_int_second[_qp](1,1) + _piezostrictive_tensor[_qp](_component,1,2) * _potential_int_second[_qp](1,2) + _piezostrictive_tensor[_qp](_component,2,0) * _potential_int_second[_qp](2,0) + _piezostrictive_tensor[_qp](_component,2,1) * _potential_int_second[_qp](2,1) + _piezostrictive_tensor[_qp](_component,2,2) * _potential_int_second[_qp](2,2));
}


Real
ConversePiezoelectricStrain::computeQpJacobian()
{
  return _test[_i][_qp] * (_piezostrictive_tensor[_qp](_component,0,0) * _second_phi[_j][_qp](0,0) + _piezostrictive_tensor[_qp](_component,0,1) * _second_phi[_j][_qp](0,1) + _piezostrictive_tensor[_qp](_component,0,2) * _second_phi[_j][_qp](0,2) + _piezostrictive_tensor[_qp](_component,1,0) * _second_phi[_j][_qp](1,0) + _piezostrictive_tensor[_qp](_component,1,1) * _second_phi[_j][_qp](1,1) + _piezostrictive_tensor[_qp](_component,1,2) * _second_phi[_j][_qp](1,2) + _piezostrictive_tensor[_qp](_component,2,0) * _second_phi[_j][_qp](2,0) + _piezostrictive_tensor[_qp](_component,2,1) * _second_phi[_j][_qp](2,1) + _piezostrictive_tensor[_qp](_component,2,2) * _second_phi[_j][_qp](2,2));
}

