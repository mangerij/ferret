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

#include "ConversePiezoelectricStrain.h"
#include "ComputePiezostrictiveTensor.h"

class ConversePiezoelectricStrain;

registerMooseObject("FerretApp", ConversePiezoelectricStrain);

InputParameters ConversePiezoelectricStrain::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates the residual for additional piezoelectric strain arising in the conditions for mechanical equilibrium.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("potential_E_int", "The electrostatic potential");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

ConversePiezoelectricStrain::ConversePiezoelectricStrain(const InputParameters & parameters)
  :Kernel(parameters),
   _piezostrictive_tensor(getMaterialProperty<RankThreeTensor>("piezostrictive_tensor")),
   _component(getParam<unsigned int>("component")),
   _potential_E_int_var(coupled("potential_E_int")),
   _potential_E_int(coupledValue("potential_E_int")),
   _potential_E_int_grad(coupledGradient("potential_E_int")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
ConversePiezoelectricStrain::computeQpResidual()
{
  Real sum = 0.0;
  for (unsigned int j = 0; j < 3; ++j)
  {
    sum += _grad_test[_i][_qp](j) * ((_potential_E_int_grad[_qp](0) * _piezostrictive_tensor[_qp](0,j,_component)) + (_potential_E_int_grad[_qp](1) * _piezostrictive_tensor[_qp](1,j,_component)) + (_potential_E_int_grad[_qp](2) * _piezostrictive_tensor[_qp](2,j,_component)));
  }
  return sum;
}

Real
ConversePiezoelectricStrain::computeQpJacobian()
{
  return 0.0;
}

Real
ConversePiezoelectricStrain::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real sum = 0.0;
  if (jvar == _potential_E_int_var)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      sum += _grad_test[_i][_qp](j) * ((_grad_phi[_j][_qp](0) * _piezostrictive_tensor[_qp](0,j,_component)) + (_grad_phi[_j][_qp](1) * _piezostrictive_tensor[_qp](1,j,_component)) + (_grad_phi[_j][_qp](2) * _piezostrictive_tensor[_qp](2,j,_component)));
    }
    return sum;
  }
  else
  {
    return 0.0;
  }
}
