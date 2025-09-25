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

#include "ElectrostrictiveCouplingPolarDerivativeRefactor.h"

class ElectrostrictiveCouplingPolarDerivativeRefactor;

registerMooseObject("FerretApp", ElectrostrictiveCouplingPolarDerivativeRefactor);

InputParameters ElectrostrictiveCouplingPolarDerivativeRefactor::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to the variation w.r.t polarization of the electrostrictive coupling energy. Note: for cubic parent phase only.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  return params;
}

ElectrostrictiveCouplingPolarDerivativeRefactor::ElectrostrictiveCouplingPolarDerivativeRefactor(const InputParameters & parameters)
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
   _Q44(getMaterialProperty<Real>("Q44")),
   _strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_strain"))
{
}

Real
ElectrostrictiveCouplingPolarDerivativeRefactor::computeQpResidual()
{
    float q11 = _C11[_qp]*_Q11[_qp] + 2*_C12[_qp]*_Q12[_qp];
    float q12 = _C11[_qp]*_Q12[_qp] + _C12[_qp]*(_Q11[_qp] + _Q12[_qp]);
    float q44 = 2*_C44[_qp]*_Q44[_qp];
  if (_component == 0)
  {
      return _test[_i][_qp]*(2*_polar_x[_qp]*((q11*_strain[_qp](0, 0) + q12*_strain[_qp](1, 1) + q12*_strain[_qp](2, 2)))
                             + 2*q44*(_polar_y[_qp]*_strain[_qp](0, 1) + _polar_z[_qp]*_strain[_qp](0, 2)));
  }
  else if (_component == 1)
  {
      return _test[_i][_qp]*(2*_polar_y[_qp]*((q11*_strain[_qp](1, 1) + q12*_strain[_qp](0, 0) + q12*_strain[_qp](2, 2)))
                             + 2*q44*(_polar_x[_qp]*_strain[_qp](0, 1) + _polar_z[_qp]*_strain[_qp](1, 2)));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp]*(2*_polar_z[_qp]*((q11*_strain[_qp](2, 2) + q12*_strain[_qp](0, 0) + q12*_strain[_qp](1, 1)))
                           + 2*q44*(_polar_y[_qp]*_strain[_qp](1, 2) + _polar_x[_qp]*_strain[_qp](0, 2)));
  }
  else
    return 0.0;
}
//
//Real
//ElectrostrictiveCouplingPolarDerivativeRefactor::computeQpJacobian()
//{
//  if (_component == 0)
//  {
//    return (-2*_phi[_j][_qp]*(_C12[_qp]*(-(_u_y_grad[_qp](1)*_Q11[_qp]) - _u_z_grad[_qp](2)*_Q11[_qp] - 2*_u_x_grad[_qp](0)*_Q12[_qp] - _u_y_grad[_qp](1)*_Q12[_qp] - _u_z_grad[_qp](2)*_Q12[_qp] + 12*Utility::pow<2>(_polar_x[_qp])*_Q11[_qp]*_Q12[_qp] + 6*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_Q12[_qp]) +
//        Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp])) + Utility::pow<2>(_polar_z[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp]))) +
//     _C11[_qp]*(-(_u_x_grad[_qp](0)*_Q11[_qp]) + 3*Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*Utility::pow<2>(_Q12[_qp])) + _Q12[_qp]*(-_u_y_grad[_qp](1) - _u_z_grad[_qp](2) + (Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*(2*_Q11[_qp] + _Q12[_qp]))) +
//     8*_C44[_qp]*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*Utility::pow<2>(_Q44[_qp]))*_test[_i][_qp]);
//  }
//  else if (_component == 1)
//  {
//    return (-2*_phi[_j][_qp]*(_C12[_qp]*(-(_u_x_grad[_qp](0)*_Q11[_qp]) - _u_z_grad[_qp](2)*_Q11[_qp] - _u_x_grad[_qp](0)*_Q12[_qp] - 2*_u_y_grad[_qp](1)*_Q12[_qp] - _u_z_grad[_qp](2)*_Q12[_qp] + 12*Utility::pow<2>(_polar_y[_qp])*_Q11[_qp]*_Q12[_qp] + 6*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_Q12[_qp]) +
//        Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp])) + Utility::pow<2>(_polar_z[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp]))) +
//     _C11[_qp]*(-(_u_y_grad[_qp](1)*_Q11[_qp]) + 3*Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*Utility::pow<2>(_Q12[_qp])) + _Q12[_qp]*(-_u_x_grad[_qp](0) - _u_z_grad[_qp](2) + (Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*(2*_Q11[_qp] + _Q12[_qp]))) +
//     8*_C44[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp]))*Utility::pow<2>(_Q44[_qp]))*_test[_i][_qp]);
//  }
//  else if (_component == 2)
//  {
//    return (-2*_phi[_j][_qp]*(_C12[_qp]*(-(_u_x_grad[_qp](0)*_Q11[_qp]) - _u_y_grad[_qp](1)*_Q11[_qp] - _u_x_grad[_qp](0)*_Q12[_qp] - _u_y_grad[_qp](1)*_Q12[_qp] - 2*_u_z_grad[_qp](2)*_Q12[_qp] + 12*Utility::pow<2>(_polar_z[_qp])*_Q11[_qp]*_Q12[_qp] + 6*Utility::pow<2>(_polar_z[_qp])*Utility::pow<2>(_Q12[_qp]) +
//        Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp])) + Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp]))) +
//     _C11[_qp]*(-(_u_z_grad[_qp](2)*_Q11[_qp]) + 3*Utility::pow<2>(_polar_z[_qp])*(Utility::pow<2>(_Q11[_qp]) + 2*Utility::pow<2>(_Q12[_qp])) + _Q12[_qp]*(-_u_x_grad[_qp](0) - _u_y_grad[_qp](1) + (Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*(2*_Q11[_qp] + _Q12[_qp]))) +
//     8*_C44[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<2>(_Q44[_qp]))*_test[_i][_qp]);
//  }
//  else
//    return 0.0;
//}
//
//Real
//ElectrostrictiveCouplingPolarDerivativeRefactor::computeQpOffDiagJacobian(unsigned int jvar)
//{
//  if (_component == 0)
//  {
//    if (jvar == _polar_y_var)
//    {
//      return (-4*_phi[_j][_qp]*(_C11[_qp]*_polar_x[_qp]*_polar_y[_qp]*_Q12[_qp]*(2*_Q11[_qp] + _Q12[_qp]) + _C12[_qp]*_polar_x[_qp]*_polar_y[_qp]*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp])) - _C44[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0))*_Q44[_qp] + 8*_C44[_qp]*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_Q44[_qp]))*_test[_i][_qp]);
//    }
//    else if (jvar == _polar_z_var)
//    {
//      return (-4*_phi[_j][_qp]*(_C11[_qp]*_polar_x[_qp]*_polar_z[_qp]*_Q12[_qp]*(2*_Q11[_qp] + _Q12[_qp]) + _C12[_qp]*_polar_x[_qp]*_polar_z[_qp]*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp])) - _C44[_qp]*(_u_x_grad[_qp](2) + _u_z_grad[_qp](0))*_Q44[_qp] + 8*_C44[_qp]*_polar_x[_qp]*_polar_z[_qp]*Utility::pow<2>(_Q44[_qp]))*_test[_i][_qp]);
//    }
//    else if (jvar == _u_x_var)
//    {
//      return (2*(_C11[_qp]*_grad_phi[_j][_qp](0)*_polar_x[_qp]*_Q11[_qp] + 2*_C12[_qp]*_grad_phi[_j][_qp](0)*_polar_x[_qp]*_Q12[_qp] + 2*_C44[_qp]*(_grad_phi[_j][_qp](1)*_polar_y[_qp] + _grad_phi[_j][_qp](2)*_polar_z[_qp])*_Q44[_qp])*_test[_i][_qp]);
//    }
//    else if (jvar == _u_y_var)
//    {
//      return (-2*(-(_C11[_qp]*_grad_phi[_j][_qp](1)*_polar_x[_qp]*_Q12[_qp]) + _C12[_qp]*_polar_x[_qp]*(-(_grad_phi[_j][_qp](1)*_Q11[_qp]) - _grad_phi[_j][_qp](1)*_Q12[_qp]) - 2*_C44[_qp]*_grad_phi[_j][_qp](0)*_polar_y[_qp]*_Q44[_qp])*_test[_i][_qp]);
//    }
//    else if (jvar == _u_z_var)
//    {
//      return (-2*(-(_C11[_qp]*_grad_phi[_j][_qp](2)*_polar_x[_qp]*_Q12[_qp]) - _C12[_qp]*_grad_phi[_j][_qp](2)*_polar_x[_qp]*(_Q11[_qp] + _Q12[_qp]) - 2*_C44[_qp]*_grad_phi[_j][_qp](0)*_polar_z[_qp]*_Q44[_qp])*_test[_i][_qp]);
//    }
//    else
//    {
//      return 0.0;
//    }
//  }
//  else if (_component == 1)
//  {
//    if (jvar == _polar_x_var)
//    {
//      return (-4*_phi[_j][_qp]*(_C11[_qp]*_polar_x[_qp]*_polar_y[_qp]*_Q12[_qp]*(2*_Q11[_qp] + _Q12[_qp]) + _C12[_qp]*_polar_x[_qp]*_polar_y[_qp]*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp])) - _C44[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0))*_Q44[_qp] + 8*_C44[_qp]*_polar_x[_qp]*_polar_y[_qp]*Utility::pow<2>(_Q44[_qp]))*_test[_i][_qp]);
//    }
//    else if (jvar == _polar_z_var)
//    {
//      return (-4*_phi[_j][_qp]*(_C11[_qp]*_polar_y[_qp]*_polar_z[_qp]*_Q12[_qp]*(2*_Q11[_qp] + _Q12[_qp]) + _C12[_qp]*_polar_y[_qp]*_polar_z[_qp]*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp])) - _C44[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))*_Q44[_qp] + 8*_C44[_qp]*_polar_y[_qp]*_polar_z[_qp]*Utility::pow<2>(_Q44[_qp]))*_test[_i][_qp]);
//    }
//    else if (jvar == _u_x_var)
//    {
//      return (2*(_C11[_qp]*_grad_phi[_j][_qp](0)*_polar_y[_qp]*_Q12[_qp] + _C12[_qp]*_grad_phi[_j][_qp](0)*_polar_y[_qp]*(_Q11[_qp] + _Q12[_qp]) + 2*_C44[_qp]*_grad_phi[_j][_qp](1)*_polar_x[_qp]*_Q44[_qp])*_test[_i][_qp]);
//    }
//    else if (jvar == _u_y_var)
//    {
//      return (2*(_C11[_qp]*_grad_phi[_j][_qp](1)*_polar_y[_qp]*_Q11[_qp] + 2*_C12[_qp]*_grad_phi[_j][_qp](1)*_polar_y[_qp]*_Q12[_qp] + 2*_C44[_qp]*(_grad_phi[_j][_qp](0)*_polar_x[_qp] + _grad_phi[_j][_qp](2)*_polar_z[_qp])*_Q44[_qp])*_test[_i][_qp]);
//    }
//    else if (jvar == _u_z_var)
//    {
//      return (2*(_C11[_qp]*_grad_phi[_j][_qp](2)*_polar_y[_qp]*_Q12[_qp] + _C12[_qp]*_grad_phi[_j][_qp](2)*_polar_y[_qp]*(_Q11[_qp] + _Q12[_qp]) + 2*_C44[_qp]*_grad_phi[_j][_qp](1)*_polar_z[_qp]*_Q44[_qp])*_test[_i][_qp]);
//    }
//    else
//    {
//      return 0.0;
//    }
//  }
//  else if (_component == 2)
//  {
//    if (jvar == _polar_x_var)
//    {
//      return (-4*_phi[_j][_qp]*(_C11[_qp]*_polar_x[_qp]*_polar_z[_qp]*_Q12[_qp]*(2*_Q11[_qp] + _Q12[_qp]) + _C12[_qp]*_polar_x[_qp]*_polar_z[_qp]*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp])) - _C44[_qp]*(_u_x_grad[_qp](2) + _u_z_grad[_qp](0))*_Q44[_qp] + 8*_C44[_qp]*_polar_x[_qp]*_polar_z[_qp]*Utility::pow<2>(_Q44[_qp]))*_test[_i][_qp]);
//    }
//    else if (jvar == _polar_y_var)
//    {
//      return (-4*_phi[_j][_qp]*(_C11[_qp]*_polar_y[_qp]*_polar_z[_qp]*_Q12[_qp]*(2*_Q11[_qp] + _Q12[_qp]) + _C12[_qp]*_polar_y[_qp]*_polar_z[_qp]*(Utility::pow<2>(_Q11[_qp]) + 2*_Q11[_qp]*_Q12[_qp] + 3*Utility::pow<2>(_Q12[_qp])) - _C44[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))*_Q44[_qp] + 8*_C44[_qp]*_polar_y[_qp]*_polar_z[_qp]*Utility::pow<2>(_Q44[_qp]))*_test[_i][_qp]);
//    }
//    else if (jvar == _u_x_var)
//    {
//      return (2*(_C11[_qp]*_grad_phi[_j][_qp](0)*_polar_z[_qp]*_Q12[_qp] + _C12[_qp]*_grad_phi[_j][_qp](0)*_polar_z[_qp]*(_Q11[_qp] + _Q12[_qp]) + 2*_C44[_qp]*_grad_phi[_j][_qp](2)*_polar_x[_qp]*_Q44[_qp])*_test[_i][_qp]);
//    }
//    else if (jvar == _u_y_var)
//    {
//      return (2*(_C11[_qp]*_grad_phi[_j][_qp](1)*_polar_z[_qp]*_Q12[_qp] + _C12[_qp]*_grad_phi[_j][_qp](1)*_polar_z[_qp]*(_Q11[_qp] + _Q12[_qp]) + 2*_C44[_qp]*_grad_phi[_j][_qp](2)*_polar_y[_qp]*_Q44[_qp])*_test[_i][_qp]);
//    }
//    else if (jvar == _u_z_var)
//    {
//      return (2*(_C11[_qp]*_grad_phi[_j][_qp](2)*_polar_z[_qp]*_Q11[_qp] + 2*_C12[_qp]*_grad_phi[_j][_qp](2)*_polar_z[_qp]*_Q12[_qp] + 2*_C44[_qp]*(_grad_phi[_j][_qp](0)*_polar_x[_qp] + _grad_phi[_j][_qp](1)*_polar_y[_qp])*_Q44[_qp])*_test[_i][_qp]);
//    }
//    else
//    {
//      return 0.0;
//    }
//  }
//  else
//    return 0.0;
//}
