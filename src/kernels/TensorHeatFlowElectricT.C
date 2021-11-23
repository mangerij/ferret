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

//Code by Dharma Raj Basaula 2021

#include "TensorHeatFlowElectricT.h"
#include "Material.h"

registerMooseObject("FerretApp", TensorHeatFlowElectricT);

InputParameters TensorHeatFlowElectricT::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to modified ohm's law");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("potential_E_int", "electrical potential");
  params.addRequiredCoupledVar("T", "temperature");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

TensorHeatFlowElectricT::TensorHeatFlowElectricT(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _potential_E_int_var(coupled("potential_E_int")),
   _potential_E_int(coupledValue("potential_E_int")),
   _potential_E_int_grad(coupledGradient("potential_E_int")),
   _T_var(coupled("T")),
   _T(coupledValue("T")),
   _T_grad(coupledGradient("T")),
   _thC_tensor(getMaterialProperty<RankTwoTensor>("thC_tensor")), //added for tensor application
   _ecC_tensor(getMaterialProperty<RankTwoTensor>("ecC_tensor")),
   _sbC_tensor(getMaterialProperty<RankTwoTensor>("sbC_tensor")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
TensorHeatFlowElectricT::computeQpResidual()
//include tensor properties//
{
  Real sum = 0.0;
  for (unsigned int i = 0, k = 0, l = 0; i < 3 && k < 3 && l < 3; ++i, ++k, ++l)
  {
    sum += ((-_grad_test[_i][_qp](0)) * (-_thC_tensor[_qp](i,0) * _T_grad[_qp](0)) +
          (-_grad_test[_i][_qp](0)) * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,l) * _sbC_tensor[_qp](l,0) * _T[_qp] * _T_grad[_qp](0)) +
          (-_grad_test[_i][_qp](0)) * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,0) * _T[_qp] * _potential_E_int_grad[_qp](0)) +
          _test[_i][_qp] * (-_ecC_tensor[_qp](i,0) * _potential_E_int_grad[_qp](0) * _potential_E_int_grad[_qp](0)) +
          _test[_i][_qp] * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,0) * _T_grad[_qp](0) * _potential_E_int_grad[_qp](0))) * _len_scale +
         ((-_grad_test[_i][_qp](1)) * (-_thC_tensor[_qp](i,1) * _T_grad[_qp](1)) +
          (-_grad_test[_i][_qp](1)) * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,l) * _sbC_tensor[_qp](l,1) * _T[_qp] * _T_grad[_qp](1)) +
          (-_grad_test[_i][_qp](1)) * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,1) * _T[_qp] * _potential_E_int_grad[_qp](1)) +
          _test[_i][_qp] * (-_ecC_tensor[_qp](i,1) * _potential_E_int_grad[_qp](1) * _potential_E_int_grad[_qp](1)) +
          _test[_i][_qp] * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,1) * _T_grad[_qp](1) * _potential_E_int_grad[_qp](1))) * _len_scale +
         ((-_grad_test[_i][_qp](2)) * (-_thC_tensor[_qp](i,2) * _T_grad[_qp](2)) +
          (-_grad_test[_i][_qp](2)) * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,l) * _sbC_tensor[_qp](l,2) * _T[_qp] * _T_grad[_qp](2)) +
          (-_grad_test[_i][_qp](2)) * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,2) * _T[_qp] * _potential_E_int_grad[_qp](2)) +
          _test[_i][_qp] * (-_ecC_tensor[_qp](i,2) * _potential_E_int_grad[_qp](2) * _potential_E_int_grad[_qp](2)) +
          _test[_i][_qp] * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,2) * _T_grad[_qp](2) * _potential_E_int_grad[_qp](2))) * _len_scale;
}
return sum;
}

Real
TensorHeatFlowElectricT::computeQpJacobian()
 //   include tensor properties//
      {
        Real sum = 0.0;
        for (unsigned int i = 0, k = 0, l = 0; i < 3 && k <3 && l < 3; ++i, ++k, ++l)
        {
         sum += ((-_grad_test[_i][_qp](0) * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,0) * _phi[_j][_qp] * _potential_E_int_grad[_qp](0)) +
                (-_grad_test[_i][_qp](0)) * (- _thC_tensor[_qp](i,0) * _grad_phi[_j][_qp](0)) +
                (-_grad_test[_i][_qp](0)) * (-_sbC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,l) * _ecC_tensor[_qp](l,0) * (_phi[_j][_qp] * _T_grad[_qp](0) + _T[_qp] * _grad_phi[_j][_qp](0))) +
                _test[_i][_qp] * (-_sbC_tensor[_qp](i,k) * _ecC_tensor[_qp](k,0) * _grad_phi[_j][_qp](0) * _potential_E_int_grad[_qp](0))) +
               ((-_grad_test[_i][_qp](1)) * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,1) * _phi[_j][_qp] * _potential_E_int_grad[_qp](1)) +
                (-_grad_test[_i][_qp](1)) * (- _thC_tensor[_qp](i,1) * _grad_phi[_j][_qp](1)) +
                (-_grad_test[_i][_qp](1)) * (-_sbC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,l) * _ecC_tensor[_qp](l,1) * (_phi[_j][_qp] * _T_grad[_qp](1) + _T[_qp] * _grad_phi[_j][_qp](1))) +
                _test[_i][_qp] * (-_sbC_tensor[_qp](i,k) * _ecC_tensor[_qp](k,1) * _grad_phi[_j][_qp](1) * _potential_E_int_grad[_qp](1))) +
               ((-_grad_test[_i][_qp](2)) * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,2) * _phi[_j][_qp] * _potential_E_int_grad[_qp](2)) +
                (-_grad_test[_i][_qp](2)) * (- _thC_tensor[_qp](i,2) * _grad_phi[_j][_qp](2)) +
                (-_grad_test[_i][_qp](2)) * (-_sbC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,l) * _ecC_tensor[_qp](l,2) * (_phi[_j][_qp] * _T_grad[_qp](2) + _T[_qp] * _grad_phi[_j][_qp](2))) +
                _test[_i][_qp] * (-_sbC_tensor[_qp](i,k) * _ecC_tensor[_qp](k,2) * _grad_phi[_j][_qp](2) * _potential_E_int_grad[_qp](2)))) * _len_scale;

            }
            return sum;
        }

Real
TensorHeatFlowElectricT::computeQpOffDiagJacobian(unsigned int jvar)
 // include tensor properties//
{
  if(jvar == _potential_E_int_var)
  {
    Real sum = 0.0;
    for (unsigned int i = 0, k = 0; i < 3 && k < 3; ++i, ++k)
   {
    sum += ((-_grad_test[_i][_qp](0)) * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,0) * _T[_qp] * _grad_phi[_j][_qp](0)) -
               _test[_i][_qp] * (2.0 * _ecC_tensor[_qp](i,0) * _potential_E_int_grad[_qp](0) * _grad_phi[_j][_qp](0) + _ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,0) * _T_grad[_qp](0) *_grad_phi[_j][_qp](0)) +
               (-_grad_test[_i][_qp](1)) * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,1) * _T[_qp] * _grad_phi[_j][_qp](1)) -
               _test[_i][_qp] * (2.0 * _ecC_tensor[_qp](i,1) * _potential_E_int_grad[_qp](1) * _grad_phi[_j][_qp](1) + _ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,1) * _T_grad[_qp](1) *_grad_phi[_j][_qp](1))  +
               (-_grad_test[_i][_qp](2)) * (-_ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,2) * _T[_qp] * _grad_phi[_j][_qp](2)) -
               _test[_i][_qp] * (2.0 * _ecC_tensor[_qp](i,2) * _potential_E_int_grad[_qp](2) * _grad_phi[_j][_qp](2) + _ecC_tensor[_qp](i,k) * _sbC_tensor[_qp](k,2) * _T_grad[_qp](2) *_grad_phi[_j][_qp](2))) * _len_scale;

    }
    return sum;
  }
    else
    {
      return 0.0;
    }
  }
