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

#include "RotostrictiveCouplingDispDerivative.h"
#include "libmesh/utility.h"

class RotostrictiveCouplingDispDerivative;

registerMooseObject("FerretApp", RotostrictiveCouplingDispDerivative);

InputParameters RotostrictiveCouplingDispDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the afd vector field");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the afd vector field");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the afd vector field");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

RotostrictiveCouplingDispDerivative::RotostrictiveCouplingDispDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _antiferrodis_A_x_var(coupled("antiferrodis_A_x")),
   _antiferrodis_A_y_var(coupled("antiferrodis_A_y")),
   _antiferrodis_A_z_var(coupled("antiferrodis_A_z")),
   _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
   _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
   _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
   _C11(getMaterialProperty<Real>("C11")),
   _C12(getMaterialProperty<Real>("C12")),
   _C44(getMaterialProperty<Real>("C44")),
   _R11(getMaterialProperty<Real>("R11")),
   _R12(getMaterialProperty<Real>("R12")),
   _R44(getMaterialProperty<Real>("R44"))
{
}

Real
RotostrictiveCouplingDispDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return ((_C11[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R11[_qp] + _C12[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R11[_qp] + _C12[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R11[_qp] + 
      2*_C12[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R12[_qp] + _C11[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R12[_qp] + _C12[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R12[_qp] + 
      _C11[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R12[_qp] + _C12[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R12[_qp])*_grad_test[_i][_qp](0) + 4*_C44[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_R44[_qp]*_grad_test[_i][_qp](1) + 
   4*_C44[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_R44[_qp]*_grad_test[_i][_qp](2));
  }
  else if (_component == 1)
  {
    return (4*_C44[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_R44[_qp]*_grad_test[_i][_qp](0) + (_C12[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R11[_qp] + _C11[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R11[_qp] + 
      _C12[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R11[_qp] + _C11[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R12[_qp] + _C12[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R12[_qp] + 
      2*_C12[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R12[_qp] + _C11[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R12[_qp] + _C12[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R12[_qp])*_grad_test[_i][_qp](1) + 
   4*_C44[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_R44[_qp]*_grad_test[_i][_qp](2));
  }
  else if (_component == 2)
  {
    return (4*_C44[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_R44[_qp]*_grad_test[_i][_qp](0) + 4*_C44[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_R44[_qp]*_grad_test[_i][_qp](1) + 
   (_C12[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R11[_qp] + _C12[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R11[_qp] + _C11[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R11[_qp] + 
      _C11[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R12[_qp] + _C12[_qp]*Utility::pow<2>(_antiferrodis_A_x[_qp])*_R12[_qp] + _C11[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R12[_qp] + 
      _C12[_qp]*Utility::pow<2>(_antiferrodis_A_y[_qp])*_R12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_antiferrodis_A_z[_qp])*_R12[_qp])*_grad_test[_i][_qp](2));
  }
  else
    return 0.0;
}

Real
RotostrictiveCouplingDispDerivative::computeQpJacobian()
{
  return 0.0;
}

Real
RotostrictiveCouplingDispDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return _phi[_j][_qp] * (2*_C11[_qp]*_antiferrodis_A_x[_qp]*_R11[_qp]*_grad_test[_i][_qp](0) + 4*_C12[_qp]*_antiferrodis_A_x[_qp]*_R12[_qp]*_grad_test[_i][_qp](0) + 4*_C44[_qp]*_antiferrodis_A_y[_qp]*_R44[_qp]*_grad_test[_i][_qp](1) + 4*_C44[_qp]*_antiferrodis_A_z[_qp]*_R44[_qp]*_grad_test[_i][_qp](2));
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _phi[_j][_qp] * (2*_C12[_qp]*_antiferrodis_A_y[_qp]*_R11[_qp]*_grad_test[_i][_qp](0) + 2*_C11[_qp]*_antiferrodis_A_y[_qp]*_R12[_qp]*_grad_test[_i][_qp](0) + 2*_C12[_qp]*_antiferrodis_A_y[_qp]*_R12[_qp]*_grad_test[_i][_qp](0) + 4*_C44[_qp]*_antiferrodis_A_x[_qp]*_R44[_qp]*_grad_test[_i][_qp](1));
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _phi[_j][_qp] * (2*_C12[_qp]*_antiferrodis_A_z[_qp]*_R11[_qp]*_grad_test[_i][_qp](0) + 2*_C11[_qp]*_antiferrodis_A_z[_qp]*_R12[_qp]*_grad_test[_i][_qp](0) + 2*_C12[_qp]*_antiferrodis_A_z[_qp]*_R12[_qp]*_grad_test[_i][_qp](0) + 4*_C44[_qp]*_antiferrodis_A_x[_qp]*_R44[_qp]*_grad_test[_i][_qp](2));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return _phi[_j][_qp] * (4*_C44[_qp]*_antiferrodis_A_y[_qp]*_R44[_qp]*_grad_test[_i][_qp](0) + 2*_C12[_qp]*_antiferrodis_A_x[_qp]*_R11[_qp]*_grad_test[_i][_qp](1) + 2*_C11[_qp]*_antiferrodis_A_x[_qp]*_R12[_qp]*_grad_test[_i][_qp](1) + 2*_C12[_qp]*_antiferrodis_A_x[_qp]*_R12[_qp]*_grad_test[_i][_qp](1));
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _phi[_j][_qp] * (4*_C44[_qp]*_antiferrodis_A_x[_qp]*_R44[_qp]*_grad_test[_i][_qp](0) + 2*_C11[_qp]*_antiferrodis_A_y[_qp]*_R11[_qp]*_grad_test[_i][_qp](1) + 4*_C12[_qp]*_antiferrodis_A_y[_qp]*_R12[_qp]*_grad_test[_i][_qp](1) + 4*_C44[_qp]*_antiferrodis_A_z[_qp]*_R44[_qp]*_grad_test[_i][_qp](2));
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _phi[_j][_qp] * (2*_C12[_qp]*_antiferrodis_A_z[_qp]*_R11[_qp]*_grad_test[_i][_qp](1) + 2*_C11[_qp]*_antiferrodis_A_z[_qp]*_R12[_qp]*_grad_test[_i][_qp](1) + 2*_C12[_qp]*_antiferrodis_A_z[_qp]*_R12[_qp]*_grad_test[_i][_qp](1) + 4*_C44[_qp]*_antiferrodis_A_y[_qp]*_R44[_qp]*_grad_test[_i][_qp](2));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _antiferrodis_A_x_var)
    {
      return _phi[_j][_qp] * (4*_C44[_qp]*_antiferrodis_A_z[_qp]*_R44[_qp]*_grad_test[_i][_qp](0) + 2*_C12[_qp]*_antiferrodis_A_x[_qp]*_R11[_qp]*_grad_test[_i][_qp](2) + 2*_C11[_qp]*_antiferrodis_A_x[_qp]*_R12[_qp]*_grad_test[_i][_qp](2) + 2*_C12[_qp]*_antiferrodis_A_x[_qp]*_R12[_qp]*_grad_test[_i][_qp](2));
    }
    else if (jvar == _antiferrodis_A_y_var)
    {
      return _phi[_j][_qp] * (4*_C44[_qp]*_antiferrodis_A_z[_qp]*_R44[_qp]*_grad_test[_i][_qp](1) + 2*_C12[_qp]*_antiferrodis_A_y[_qp]*_R11[_qp]*_grad_test[_i][_qp](2) + 2*_C11[_qp]*_antiferrodis_A_y[_qp]*_R12[_qp]*_grad_test[_i][_qp](2) + 2*_C12[_qp]*_antiferrodis_A_y[_qp]*_R12[_qp]*_grad_test[_i][_qp](2));
    }
    else if (jvar == _antiferrodis_A_z_var)
    {
      return _phi[_j][_qp] * (4*_C44[_qp]*_antiferrodis_A_x[_qp]*_R44[_qp]*_grad_test[_i][_qp](0) + 4*_C44[_qp]*_antiferrodis_A_y[_qp]*_R44[_qp]*_grad_test[_i][_qp](1) + 2*_C11[_qp]*_antiferrodis_A_z[_qp]*_R11[_qp]*_grad_test[_i][_qp](2) + 4*_C12[_qp]*_antiferrodis_A_z[_qp]*_R12[_qp]*_grad_test[_i][_qp](2));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
