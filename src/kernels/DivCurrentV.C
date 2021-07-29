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

#include "DivCurrentV.h"
// #include "HeatConduction.h"
#include "Material.h"

registerMooseObject("FerretApp", DivCurrentV);

template <>
InputParameters
validParams<DivCurrentV>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to modified ohm's law");
  params.addRequiredCoupledVar("potential_E_int", "electrical potential");
  params.addRequiredCoupledVar("T", "temperature");
  return params;
}

DivCurrentV::DivCurrentV(const InputParameters & parameters)
  : Kernel(parameters),
    _potential_E_int_var(coupled("potential_E_int")),
    _potential_E_int(coupledValue("potential_E_int")),
    _potential_E_int_grad(coupledGradient("potential_E_int")),
    _T_var(coupled("T")),
    _T(coupledValue("T")),
    _T_grad(coupledGradient("T")),
    _ecC(getMaterialProperty<Real>("ecC")),
    _sbC(getMaterialProperty<Real>("sbC"))
{
}

Real
DivCurrentV::computeQpResidual()

{
  return (((-_grad_test[_i][_qp](0)) * (-_ecC[_qp] * _sbC[_qp] * _T_grad[_qp](0)) +
           (-_grad_test[_i][_qp](0)) * (-_ecC[_qp] * _potential_E_int_grad[_qp](0))) +
          (-_grad_test[_i][_qp](1) * (-_ecC[_qp] * _sbC[_qp] * _T_grad[_qp](1)) +
           (-_grad_test[_i][_qp](1)) * (-_ecC[_qp] * _potential_E_int_grad[_qp](1))) +
          (-_grad_test[_i][_qp](2) * (-_ecC[_qp] * _sbC[_qp] * _T_grad[_qp](2)) +
           (-_grad_test[_i][_qp](2)) * (-_ecC[_qp] * _potential_E_int_grad[_qp](2))));
}

Real
DivCurrentV::computeQpJacobian()

{
  return ((-_grad_test[_i][_qp](0)) * (-_ecC[_qp] * _grad_phi[_j][_qp](0)) +
          (-_grad_test[_i][_qp](1)) * (-_ecC[_qp] * _grad_phi[_j][_qp](1)) +
          (-_grad_test[_i][_qp](2)) * (-_ecC[_qp] * _grad_phi[_j][_qp](2)));
}

Real
DivCurrentV::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _T_var)
  {
    return ((-_grad_test[_i][_qp](0)) * (-_ecC[_qp] * _sbC[_qp] * _grad_phi[_j][_qp](0)) +
            (-_grad_test[_i][_qp](1)) * (-_ecC[_qp] * _sbC[_qp] * _grad_phi[_j][_qp](1)) +
            (-_grad_test[_i][_qp](2)) * (-_ecC[_qp] * _sbC[_qp] * _grad_phi[_j][_qp](2)));
  }
  else
  {
    return 0.0;
  }
}
