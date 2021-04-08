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

#include "SeebeckEffect.h"

registerMooseObject("FerretApp", SeebeckEffect);

template <>
InputParameters
validParams<SeebeckEffect>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a contribution due to nabla.j = 0");
  params.addParam<MaterialPropertyName>("sbC", "seebeck_coefficient");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction the variable "
                                        "this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("potential_E_int", "electrical potential");
  params.addRequiredCoupledVar("T", "temperature");
  // params.addParam<MaterialPropertyName>("ecC", "Electrical Conductivity", "Property name of the
  // electrical conductivity material property");
  params.addParam<MaterialPropertyName>(
      "sbC", "Seebeck coefficient", "Property name of the Seebeck coefficient material property");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

SeebeckEffect::SeebeckEffect(const InputParameters & parameters)
  : Kernel(parameters),
    _component(getParam<unsigned int>("component")),
    _potential_E_int_var(coupled("potential_E_int")),
    _potential_E_int(coupledValue("potential_E_int")),
    _potential_E_int_grad(coupledGradient("potential_E_int")),
    _T_var(coupled("T")),
    _T(coupledValue("T")),
    _T_grad(coupledGradient("T")),
    // _ecC(getMaterialProperty<Real>("ecC")),
    _sbC(getMaterialProperty<Real>("sbC")),
    _len_scale(getParam<Real>("len_scale"))
{
}

Real
SeebeckEffect::computeQpResidual()
{
  return -_grad_test[_i][_qp](_component) *
         (_potential_E_int_grad[_qp](_component) + _sbC[_qp] * _T_grad[_qp](_component)) *
         _len_scale;
}

Real
SeebeckEffect::computeQpJacobian()
{
  // return - _grad_test[_i][_qp](_component) * (_potential_E_int_grad[_qp](_component) + _sbC[_qp]
  // * _grad_phi[_j][_qp](_component)) * _len_scale;
  return -_grad_test[_i][_qp](_component) *
         (_grad_phi[_j][_qp](_component) + _sbC[_qp] * _T_grad[_qp](_component)) * _len_scale;
}

Real
SeebeckEffect::computeQpOffDiagJacobian(unsigned int jvar)
{
  // if(jvar == _potential_E_int_var)
  //  {
  //   return (-_grad_test[_i][_qp](_component)) * _grad_phi[_j][_qp](_component) * _len_scale;
  //
  //   }
  //   else
  //   {
  //     return 0.0;
  //   }
  if (jvar == _T_var)
  {
    return (-_grad_test[_i][_qp](_component)) * _sbC[_qp] * _grad_phi[_j][_qp](_component) *
           _len_scale;
  }
  else
  {
    return 0.0;
  }
}
