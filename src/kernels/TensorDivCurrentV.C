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

//Code by Dharma Raj Basaula 2021

#include "TensorDivCurrentV.h"
// #include "HeatConduction.h"
#include "Material.h"

registerMooseObject("FerretApp", TensorDivCurrentV);

InputParameters TensorDivCurrentV::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to modified ohm's law");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("potential_E_int", "electrical potential");
  params.addRequiredCoupledVar("T", "temperature");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

TensorDivCurrentV::TensorDivCurrentV(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _potential_E_int_var(coupled("potential_E_int")),
   _potential_E_int(coupledValue("potential_E_int")),
   _potential_E_int_grad(coupledGradient("potential_E_int")),
   _T_var(coupled("T")),
   _T(coupledValue("T")),
   _T_grad(coupledGradient("T")),
   _ecC_tensor(getMaterialProperty<RankTwoTensor>("ecC_tensor")), //added for tensor application
   _sbC_tensor(getMaterialProperty<RankTwoTensor>("sbC_tensor")), //added for tensor application
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
TensorDivCurrentV::computeQpResidual()
 //include tensor//

 {
   Real sum = 0.0;
   for (unsigned int j = 0, k = 0; j < 3 && k < 3; ++j, ++k)
   {
     sum += (((-_grad_test[_i][_qp](0)) * (-_ecC_tensor[_qp](j,0) * _sbC_tensor[_qp](0,k) * _T_grad[_qp](0)) +
            (-_grad_test[_i][_qp](0)) * (-_ecC_tensor[_qp](j,0) * _potential_E_int_grad[_qp](0))) +
            (-_grad_test[_i][_qp](1) * (-_ecC_tensor[_qp](j,1) * _sbC_tensor[_qp](1,k) * _T_grad[_qp](1)) +
            (-_grad_test[_i][_qp](1)) * (-_ecC_tensor[_qp](j,1) * _potential_E_int_grad[_qp](1))) +
            (-_grad_test[_i][_qp](2) * (-_ecC_tensor[_qp](j,2) * _sbC_tensor[_qp](2,k) * _T_grad[_qp](2)) +
            (-_grad_test[_i][_qp](2)) * (-_ecC_tensor[_qp](j,2) * _potential_E_int_grad[_qp](2)))) * _len_scale;
    }
    return sum;
  }

Real
TensorDivCurrentV::computeQpJacobian()

//include tensor//
{
  Real sum = 0.0;
  for (unsigned int j = 0; j < 3; ++j)
  {
    sum += ((-_grad_test[_i][_qp](0)) * (- _ecC_tensor[_qp](j,0) * _grad_phi[_j][_qp](0)) +
              (-_grad_test[_i][_qp](1)) * (-_ecC_tensor[_qp](j,1) * _grad_phi[_j][_qp](1)) +
              (-_grad_test[_i][_qp](2)) * (-_ecC_tensor[_qp](j,2) * _grad_phi[_j][_qp](2))) * _len_scale;
   }
   return sum;
 }


Real
TensorDivCurrentV::computeQpOffDiagJacobian(unsigned int jvar)
  //include tensor//

  {
    if(jvar == _T_var)
    {
      Real sum = 0.0;
      for (unsigned int j = 0, k = 0; j < 3 && k < 3; ++j, ++k)
     {
      sum += ((-_grad_test[_i][_qp](0)) * (-_ecC_tensor[_qp](j,0) * _sbC_tensor[_qp](0,k) * _grad_phi[_j][_qp](0)) +
             (-_grad_test[_i][_qp](1)) * (-_ecC_tensor[_qp](j,1) * _sbC_tensor[_qp](1,k) * _grad_phi[_j][_qp](1)) +
             (-_grad_test[_i][_qp](2)) * (-_ecC_tensor[_qp](j,2) * _sbC_tensor[_qp](2,k) * _grad_phi[_j][_qp](2))) * _len_scale;

      }
      return sum;
    }
      else
      {
        return 0.0;
      }
    }
