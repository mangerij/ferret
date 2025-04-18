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

   You should have received a co_polar_y[_qp] of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "ElectricFluxTensor.h"
#include "Material.h"
#include "RankTwoTensor.h"

registerMooseObject("FerretApp", ElectricFluxTensor);

InputParameters ElectricFluxTensor::validParams()
{
InputParameters params = AuxKernel::validParams();
params.addClassDescription("Electric flux generated");
params.addRequiredCoupledVar("T", "temperature");
params.addRequiredCoupledVar("potential_E_int", "electric potential");
params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
return params;
}

ElectricFluxTensor::ElectricFluxTensor(const InputParameters & parameters)
  : AuxKernel(parameters),
    _T(coupledValue("T")),
    _T_grad(coupledGradient("T")),
    _potential_E_int(coupledValue("potential_E_int")),
    _potential_E_int_grad(coupledGradient("potential_E_int")),
    _ecC_tensor(getMaterialProperty<RankTwoTensor>("ecC_tensor")),
    _sbC_tensor(getMaterialProperty<RankTwoTensor>("sbC_tensor")),
    _component(getParam<unsigned int>("component"))
{
}

Real
ElectricFluxTensor::computeValue()
{
  Real sum = 0.0;
  for (unsigned int i = 0, j = 0, k = 0; i < 3 && j < 3 && k < 3; ++i, ++j, ++k)
  {
     sum += -_ecC_tensor[_qp](i,j) * _potential_E_int_grad[_qp](_component) -
             _sbC_tensor[_qp](i,j) * _ecC_tensor[_qp](j,k) * _T_grad[_qp](_component);
  }
  return sum;
}
