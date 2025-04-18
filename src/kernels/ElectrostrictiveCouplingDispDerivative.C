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

#include "ElectrostrictiveCouplingDispDerivative.h"
#include "libmesh/utility.h"

class ElectrostrictiveCouplingDispDerivative;

registerMooseObject("FerretApp", ElectrostrictiveCouplingDispDerivative);

InputParameters ElectrostrictiveCouplingDispDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to the spontaneous ferroelectric strain in the condition for mechanical equilibrium.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  return params;
}

ElectrostrictiveCouplingDispDerivative::ElectrostrictiveCouplingDispDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _C11(getMaterialProperty<Real>("C11")),
   _C12(getMaterialProperty<Real>("C12")),
   _C44(getMaterialProperty<Real>("C44")),
   _Q11(getMaterialProperty<Real>("Q11")),
   _Q12(getMaterialProperty<Real>("Q12")),
   _Q44(getMaterialProperty<Real>("Q44"))
{
}

Real
ElectrostrictiveCouplingDispDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return ((_C11[_qp]*Utility::pow<2>(_polar_x[_qp])*_Q11[_qp] + _C12[_qp]*Utility::pow<2>(_polar_y[_qp])*_Q11[_qp] + _C12[_qp]*Utility::pow<2>(_polar_z[_qp])*_Q11[_qp] +
      2*_C12[_qp]*Utility::pow<2>(_polar_x[_qp])*_Q12[_qp] + _C11[_qp]*Utility::pow<2>(_polar_y[_qp])*_Q12[_qp] + _C12[_qp]*Utility::pow<2>(_polar_y[_qp])*_Q12[_qp] +
      _C11[_qp]*Utility::pow<2>(_polar_z[_qp])*_Q12[_qp] + _C12[_qp]*Utility::pow<2>(_polar_z[_qp])*_Q12[_qp])*_grad_test[_i][_qp](0) + 4*_C44[_qp]*_polar_x[_qp]*_polar_y[_qp]*_Q44[_qp]*_grad_test[_i][_qp](1) +
   4*_C44[_qp]*_polar_x[_qp]*_polar_z[_qp]*_Q44[_qp]*_grad_test[_i][_qp](2));
  }
  else if (_component == 1)
  {
    return (4*_C44[_qp]*_polar_x[_qp]*_polar_y[_qp]*_Q44[_qp]*_grad_test[_i][_qp](0) + (_C12[_qp]*Utility::pow<2>(_polar_x[_qp])*_Q11[_qp] + _C11[_qp]*Utility::pow<2>(_polar_y[_qp])*_Q11[_qp] +
      _C12[_qp]*Utility::pow<2>(_polar_z[_qp])*_Q11[_qp] + _C11[_qp]*Utility::pow<2>(_polar_x[_qp])*_Q12[_qp] + _C12[_qp]*Utility::pow<2>(_polar_x[_qp])*_Q12[_qp] +
      2*_C12[_qp]*Utility::pow<2>(_polar_y[_qp])*_Q12[_qp] + _C11[_qp]*Utility::pow<2>(_polar_z[_qp])*_Q12[_qp] + _C12[_qp]*Utility::pow<2>(_polar_z[_qp])*_Q12[_qp])*_grad_test[_i][_qp](1) +
   4*_C44[_qp]*_polar_y[_qp]*_polar_z[_qp]*_Q44[_qp]*_grad_test[_i][_qp](2));
  }
  else if (_component == 2)
  {
    return (4*_C44[_qp]*_polar_x[_qp]*_polar_z[_qp]*_Q44[_qp]*_grad_test[_i][_qp](0) + 4*_C44[_qp]*_polar_y[_qp]*_polar_z[_qp]*_Q44[_qp]*_grad_test[_i][_qp](1) +
   (_C12[_qp]*Utility::pow<2>(_polar_x[_qp])*_Q11[_qp] + _C12[_qp]*Utility::pow<2>(_polar_y[_qp])*_Q11[_qp] + _C11[_qp]*Utility::pow<2>(_polar_z[_qp])*_Q11[_qp] +
      _C11[_qp]*Utility::pow<2>(_polar_x[_qp])*_Q12[_qp] + _C12[_qp]*Utility::pow<2>(_polar_x[_qp])*_Q12[_qp] + _C11[_qp]*Utility::pow<2>(_polar_y[_qp])*_Q12[_qp] +
      _C12[_qp]*Utility::pow<2>(_polar_y[_qp])*_Q12[_qp] + 2*_C12[_qp]*Utility::pow<2>(_polar_z[_qp])*_Q12[_qp])*_grad_test[_i][_qp](2));
  }
  else
    return 0.0;
}

Real
ElectrostrictiveCouplingDispDerivative::computeQpJacobian()
{
  return 0.0;
}

Real
ElectrostrictiveCouplingDispDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_x_var)
    {
      return _phi[_j][_qp] * (2*_C11[_qp]*_polar_x[_qp]*_Q11[_qp]*_grad_test[_i][_qp](0) + 4*_C12[_qp]*_polar_x[_qp]*_Q12[_qp]*_grad_test[_i][_qp](0) + 4*_C44[_qp]*_polar_y[_qp]*_Q44[_qp]*_grad_test[_i][_qp](1) + 4*_C44[_qp]*_polar_z[_qp]*_Q44[_qp]*_grad_test[_i][_qp](2));
    }
    else if (jvar == _polar_y_var)
    {
      return _phi[_j][_qp] * (2*_C12[_qp]*_polar_y[_qp]*_Q11[_qp]*_grad_test[_i][_qp](0) + 2*_C11[_qp]*_polar_y[_qp]*_Q12[_qp]*_grad_test[_i][_qp](0) + 2*_C12[_qp]*_polar_y[_qp]*_Q12[_qp]*_grad_test[_i][_qp](0) + 4*_C44[_qp]*_polar_x[_qp]*_Q44[_qp]*_grad_test[_i][_qp](1));
    }
    else if (jvar == _polar_z_var)
    {
      return _phi[_j][_qp] * (2*_C12[_qp]*_polar_z[_qp]*_Q11[_qp]*_grad_test[_i][_qp](0) + 2*_C11[_qp]*_polar_z[_qp]*_Q12[_qp]*_grad_test[_i][_qp](0) + 2*_C12[_qp]*_polar_z[_qp]*_Q12[_qp]*_grad_test[_i][_qp](0) + 4*_C44[_qp]*_polar_x[_qp]*_Q44[_qp]*_grad_test[_i][_qp](2));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _polar_x_var)
    {
      return _phi[_j][_qp] * (4*_C44[_qp]*_polar_y[_qp]*_Q44[_qp]*_grad_test[_i][_qp](0) + 2*_C12[_qp]*_polar_x[_qp]*_Q11[_qp]*_grad_test[_i][_qp](1) + 2*_C11[_qp]*_polar_x[_qp]*_Q12[_qp]*_grad_test[_i][_qp](1) + 2*_C12[_qp]*_polar_x[_qp]*_Q12[_qp]*_grad_test[_i][_qp](1));
    }
    else if (jvar == _polar_y_var)
    {
      return _phi[_j][_qp] * (4*_C44[_qp]*_polar_x[_qp]*_Q44[_qp]*_grad_test[_i][_qp](0) + 2*_C11[_qp]*_polar_y[_qp]*_Q11[_qp]*_grad_test[_i][_qp](1) + 4*_C12[_qp]*_polar_y[_qp]*_Q12[_qp]*_grad_test[_i][_qp](1) + 4*_C44[_qp]*_polar_z[_qp]*_Q44[_qp]*_grad_test[_i][_qp](2));
    }
    else if (jvar == _polar_z_var)
    {
      return _phi[_j][_qp] * (2*_C12[_qp]*_polar_z[_qp]*_Q11[_qp]*_grad_test[_i][_qp](1) + 2*_C11[_qp]*_polar_z[_qp]*_Q12[_qp]*_grad_test[_i][_qp](1) + 2*_C12[_qp]*_polar_z[_qp]*_Q12[_qp]*_grad_test[_i][_qp](1) + 4*_C44[_qp]*_polar_y[_qp]*_Q44[_qp]*_grad_test[_i][_qp](2));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _polar_x_var)
    {
      return _phi[_j][_qp] * (4*_C44[_qp]*_polar_z[_qp]*_Q44[_qp]*_grad_test[_i][_qp](0) + 2*_C12[_qp]*_polar_x[_qp]*_Q11[_qp]*_grad_test[_i][_qp](2) + 2*_C11[_qp]*_polar_x[_qp]*_Q12[_qp]*_grad_test[_i][_qp](2) + 2*_C12[_qp]*_polar_x[_qp]*_Q12[_qp]*_grad_test[_i][_qp](2));
    }
    else if (jvar == _polar_y_var)
    {
      return _phi[_j][_qp] * (4*_C44[_qp]*_polar_z[_qp]*_Q44[_qp]*_grad_test[_i][_qp](1) + 2*_C12[_qp]*_polar_y[_qp]*_Q11[_qp]*_grad_test[_i][_qp](2) + 2*_C11[_qp]*_polar_y[_qp]*_Q12[_qp]*_grad_test[_i][_qp](2) + 2*_C12[_qp]*_polar_y[_qp]*_Q12[_qp]*_grad_test[_i][_qp](2));
    }
    else if (jvar == _polar_z_var)
    {
      return _phi[_j][_qp] * (4*_C44[_qp]*_polar_x[_qp]*_Q44[_qp]*_grad_test[_i][_qp](0) + 4*_C44[_qp]*_polar_y[_qp]*_Q44[_qp]*_grad_test[_i][_qp](1) + 2*_C11[_qp]*_polar_z[_qp]*_Q11[_qp]*_grad_test[_i][_qp](2) + 4*_C12[_qp]*_polar_z[_qp]*_Q12[_qp]*_grad_test[_i][_qp](2));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
