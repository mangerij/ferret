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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "AFMInteractionCartLL.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", AFMInteractionCartLL);

InputParameters AFMInteractionCartLL::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution for the sublattice exchange in an antiferromagnet");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredParam<unsigned int>("mag_sub", "An integer corresponding to the sublattice this Kernel acts on");
  params.addRequiredCoupledVar("potential_H_int", "The internal magnetic potential variable");
  params.addRequiredCoupledVar("mag1_x", "The x component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag1_y", "The y component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag1_z", "The z component of the constrained 1st sublattice magnetization vector");
  params.addCoupledVar("mag2_x", 0.0, "The x component of the constrained 2nd sublattice magnetization vector");
  params.addCoupledVar("mag2_y", 0.0, "The y component of the constrained 2nd sublattice magnetization vector");
  params.addCoupledVar("mag2_z", 0.0, "The z component of the constrained 2nd sublattice magnetization vector");
  return params;
}

AFMInteractionCartLL::AFMInteractionCartLL(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _mag_sub(getParam<unsigned int>("mag_sub")),
   _potential_H_int_var(coupled("potential_H_int")),
   _potential_H_int(coupledValue("potential_H_int")),
   _potential_H_int_grad(coupledGradient("potential_H_int")),
   _mag1_x_var(coupled("mag1_x")),
   _mag1_y_var(coupled("mag1_y")),
   _mag1_z_var(coupled("mag1_z")),
   _mag1_x(coupledValue("mag1_x")),
   _mag1_y(coupledValue("mag1_y")),
   _mag1_z(coupledValue("mag1_z")),
   _mag2_x_var(coupled("mag2_x")),
   _mag2_y_var(coupled("mag2_y")),
   _mag2_z_var(coupled("mag2_z")),
   _mag2_x(coupledValue("mag2_x")),
   _mag2_y(coupledValue("mag2_y")),
   _mag2_z(coupledValue("mag2_z")),
   _alpha(getMaterialProperty<Real>("alpha")),
   _g0(getMaterialProperty<Real>("g0")),
   _Ms(getMaterialProperty<Real>("Ms"))
{
}


Real
AFMInteractionCartLL::computeQpResidual()
{
  if (_mag_sub == 0)
  {
    if (_component == 0)
    {
      return (_g0[_qp]*_Ms[_qp]*(_alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag1_x[_qp]*_mag1_y[_qp] - _potential_H_int_grad[_qp](1)*_mag1_z[_qp] + _potential_H_int_grad[_qp](2)*(_mag1_y[_qp] + _alpha[_qp]*_mag1_x[_qp]*_mag1_z[_qp]) - _alpha[_qp]*_potential_H_int_grad[_qp](0)*(Utility::pow<2>(_mag1_y[_qp]) + Utility::pow<2>(_mag1_z[_qp])))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 1)
    {
      return (_g0[_qp]*_Ms[_qp]*(_potential_H_int_grad[_qp](0)*(_alpha[_qp]*_mag1_x[_qp]*_mag1_y[_qp] + _mag1_z[_qp]) + _potential_H_int_grad[_qp](2)*(-_mag1_x[_qp] + _alpha[_qp]*_mag1_y[_qp]*_mag1_z[_qp]) - _alpha[_qp]*_potential_H_int_grad[_qp](1)*(Utility::pow<2>(_mag1_x[_qp]) + Utility::pow<2>(_mag1_z[_qp])))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 2)
    {
      return (_g0[_qp]*_Ms[_qp]*(-(_potential_H_int_grad[_qp](0)*_mag1_y[_qp]) - _alpha[_qp]*_potential_H_int_grad[_qp](2)*(Utility::pow<2>(_mag1_x[_qp]) + Utility::pow<2>(_mag1_y[_qp])) + _alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag1_x[_qp]*_mag1_z[_qp] + _potential_H_int_grad[_qp](1)*(_mag1_x[_qp] + _alpha[_qp]*_mag1_y[_qp]*_mag1_z[_qp]))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else
      return 0.0;
  }
  else if (_mag_sub == 1)
  {
    if (_component == 0)
    {
      return (_g0[_qp]*_Ms[_qp]*(_alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag2_x[_qp]*_mag2_y[_qp] - _potential_H_int_grad[_qp](1)*_mag2_z[_qp] + _potential_H_int_grad[_qp](2)*(_mag2_y[_qp] + _alpha[_qp]*_mag2_x[_qp]*_mag2_z[_qp]) - _alpha[_qp]*_potential_H_int_grad[_qp](0)*(Utility::pow<2>(_mag2_y[_qp]) + Utility::pow<2>(_mag2_z[_qp])))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 1)
    {
      return (_g0[_qp]*_Ms[_qp]*(_potential_H_int_grad[_qp](0)*(_alpha[_qp]*_mag2_x[_qp]*_mag2_y[_qp] + _mag2_z[_qp]) + _potential_H_int_grad[_qp](2)*(-_mag2_x[_qp] + _alpha[_qp]*_mag2_y[_qp]*_mag2_z[_qp]) - _alpha[_qp]*_potential_H_int_grad[_qp](1)*(Utility::pow<2>(_mag2_x[_qp]) + Utility::pow<2>(_mag2_z[_qp])))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 2)
    {
      return (_g0[_qp]*_Ms[_qp]*(-(_potential_H_int_grad[_qp](0)*_mag2_y[_qp]) - _alpha[_qp]*_potential_H_int_grad[_qp](2)*(Utility::pow<2>(_mag2_x[_qp]) + Utility::pow<2>(_mag2_y[_qp])) + _alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag2_x[_qp]*_mag2_z[_qp] + _potential_H_int_grad[_qp](1)*(_mag2_x[_qp] + _alpha[_qp]*_mag2_y[_qp]*_mag2_z[_qp]))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}

Real
AFMInteractionCartLL::computeQpJacobian()
{
  if (_mag_sub == 0)
  {
    if (_component == 0)
    {
      return (_g0[_qp]*_Ms[_qp]*(_alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag1_y[_qp] + _alpha[_qp]*_potential_H_int_grad[_qp](2)*_mag1_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 1)
    {
      return (_g0[_qp]*_Ms[_qp]*(_alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag1_x[_qp] + _alpha[_qp]*_potential_H_int_grad[_qp](2)*_mag1_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 2)
    {
      return (_g0[_qp]*_Ms[_qp]*(_alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag1_x[_qp] + _alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag1_y[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else
      return 0.0;
  }
  else if (_mag_sub == 1)
  {
    if (_component == 0)
    {
      return (_g0[_qp]*_Ms[_qp]*(_alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag2_y[_qp] + _alpha[_qp]*_potential_H_int_grad[_qp](2)*_mag2_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 1)
    {
      return (_g0[_qp]*_Ms[_qp]*(_alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag2_x[_qp] + _alpha[_qp]*_potential_H_int_grad[_qp](2)*_mag2_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else if (_component == 2)
    {
      return (_g0[_qp]*_Ms[_qp]*(_alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag2_x[_qp] + _alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag2_y[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}

Real
AFMInteractionCartLL::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_mag_sub == 0)
  {
    if (_component == 0)
    {
      if (jvar == _mag1_y_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(_potential_H_int_grad[_qp](2) + _alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag1_x[_qp] - 2.0*_alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag1_y[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag1_z_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(-_potential_H_int_grad[_qp](1) + _alpha[_qp]*_potential_H_int_grad[_qp](2)*_mag1_x[_qp] - 2.0*_alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag1_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag2_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_y_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_z_var)
      {
        return 0.0;
      }
      else if (jvar == _potential_H_int_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(_alpha[_qp]*_grad_phi[_j][_qp](1)*_mag1_x[_qp]*_mag1_y[_qp] - _grad_phi[_j][_qp](1)*_mag1_z[_qp] + _grad_phi[_j][_qp](2)*(_mag1_y[_qp] + _alpha[_qp]*_mag1_x[_qp]*_mag1_z[_qp]) - _alpha[_qp]*_grad_phi[_j][_qp](0)*(Utility::pow<2>(_mag1_y[_qp]) + Utility::pow<2>(_mag1_z[_qp])))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else
        return 0.0;
    }
    else if (_component == 1)
    {
      if (jvar == _mag1_x_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(-_potential_H_int_grad[_qp](2) - 2.0*_alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag1_x[_qp] + _alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag1_y[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag1_z_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(_potential_H_int_grad[_qp](0) + _alpha[_qp]*_potential_H_int_grad[_qp](2)*_mag1_y[_qp] - 2.0*_alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag1_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag2_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_y_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_z_var)
      {
        return 0.0;
      }
      else if (jvar == _potential_H_int_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(_grad_phi[_j][_qp](0)*(_alpha[_qp]*_mag1_x[_qp]*_mag1_y[_qp] + _mag1_z[_qp]) + _grad_phi[_j][_qp](2)*(-_mag1_x[_qp] + _alpha[_qp]*_mag1_y[_qp]*_mag1_z[_qp]) - _alpha[_qp]*_grad_phi[_j][_qp](1)*(Utility::pow<2>(_mag1_x[_qp]) + Utility::pow<2>(_mag1_z[_qp])))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else
        return 0.0;
    }
    else if (_component == 2)
    {
      if (jvar == _mag1_x_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(_potential_H_int_grad[_qp](1) - 2.0*_alpha[_qp]*_potential_H_int_grad[_qp](2)*_mag1_x[_qp] + _alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag1_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag1_y_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(-_potential_H_int_grad[_qp](0) - 2.0*_alpha[_qp]*_potential_H_int_grad[_qp](2)*_mag1_y[_qp] + _alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag1_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag2_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_y_var)
      {
        return 0.0;
      }
      else if (jvar == _mag2_z_var)
      {
        return 0.0;
      }
      else if (jvar == _potential_H_int_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(-(_grad_phi[_j][_qp](0)*_mag1_y[_qp]) - _alpha[_qp]*_grad_phi[_j][_qp](2)*(Utility::pow<2>(_mag1_x[_qp]) + Utility::pow<2>(_mag1_y[_qp])) + _alpha[_qp]*_grad_phi[_j][_qp](0)*_mag1_x[_qp]*_mag1_z[_qp] + _grad_phi[_j][_qp](1)*(_mag1_x[_qp] + _alpha[_qp]*_mag1_y[_qp]*_mag1_z[_qp]))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else
        return 0.0;
    }
    else
      return 0.0;
  }
  else if (_mag_sub == 1)
  {
    if (_component == 0)
    {
      if (jvar == _mag2_y_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(_potential_H_int_grad[_qp](2) + _alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag2_x[_qp] - 2.0*_alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag2_y[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag2_z_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(-_potential_H_int_grad[_qp](1) + _alpha[_qp]*_potential_H_int_grad[_qp](2)*_mag2_x[_qp] - 2.0*_alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag2_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag1_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_y_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_z_var)
      {
        return 0.0;
      }
      else if (jvar == _potential_H_int_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(_alpha[_qp]*_grad_phi[_j][_qp](1)*_mag2_x[_qp]*_mag2_y[_qp] - _grad_phi[_j][_qp](1)*_mag2_z[_qp] + _grad_phi[_j][_qp](2)*(_mag2_y[_qp] + _alpha[_qp]*_mag2_x[_qp]*_mag2_z[_qp]) - _alpha[_qp]*_grad_phi[_j][_qp](0)*(Utility::pow<2>(_mag2_y[_qp]) + Utility::pow<2>(_mag2_z[_qp])))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else
        return 0.0;
    }
    else if (_component == 1)
    {
      if (jvar == _mag2_x_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(-_potential_H_int_grad[_qp](2) - 2.0*_alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag2_x[_qp] + _alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag2_y[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag2_z_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(_potential_H_int_grad[_qp](0) + _alpha[_qp]*_potential_H_int_grad[_qp](2)*_mag2_y[_qp] - 2.0*_alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag2_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag1_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_y_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_z_var)
      {
        return 0.0;
      }
      else if (jvar == _potential_H_int_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(_grad_phi[_j][_qp](0)*(_alpha[_qp]*_mag2_x[_qp]*_mag2_y[_qp] + _mag2_z[_qp]) + _grad_phi[_j][_qp](2)*(-_mag2_x[_qp] + _alpha[_qp]*_mag2_y[_qp]*_mag2_z[_qp]) - _alpha[_qp]*_grad_phi[_j][_qp](1)*(Utility::pow<2>(_mag2_x[_qp]) + Utility::pow<2>(_mag2_z[_qp])))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else
        return 0.0;
    }
    else if (_component == 2)
    {
      if (jvar == _mag2_x_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(_potential_H_int_grad[_qp](1) - 2.0*_alpha[_qp]*_potential_H_int_grad[_qp](2)*_mag2_x[_qp] + _alpha[_qp]*_potential_H_int_grad[_qp](0)*_mag2_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag2_y_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(-_potential_H_int_grad[_qp](0) - 2.0*_alpha[_qp]*_potential_H_int_grad[_qp](2)*_mag2_y[_qp] + _alpha[_qp]*_potential_H_int_grad[_qp](1)*_mag2_z[_qp])*_phi[_j][_qp]*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else if (jvar == _mag1_x_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_y_var)
      {
        return 0.0;
      }
      else if (jvar == _mag1_z_var)
      {
        return 0.0;
      }
      else if (jvar == _potential_H_int_var)
      {
        return (_g0[_qp]*_Ms[_qp]*(-(_grad_phi[_j][_qp](0)*_mag2_y[_qp]) - _alpha[_qp]*_grad_phi[_j][_qp](2)*(Utility::pow<2>(_mag2_x[_qp]) + Utility::pow<2>(_mag2_y[_qp])) + _alpha[_qp]*_grad_phi[_j][_qp](0)*_mag2_x[_qp]*_mag2_z[_qp] + _grad_phi[_j][_qp](1)*(_mag2_x[_qp] + _alpha[_qp]*_mag2_y[_qp]*_mag2_z[_qp]))*_test[_i][_qp])/(1.0 + Utility::pow<2>(_alpha[_qp]));
      }
      else
        return 0.0;
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}
