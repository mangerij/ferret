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

#include "Transformed111KernelOp3.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", Transformed111KernelOp3);

InputParameters Transformed111KernelOp3::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates the transformed residual for the local free energy which is an eighth order expansion in the polarization.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  params.addRequiredCoupledVar("order_param_x", "The x component of the transformed order parameter");
  params.addRequiredCoupledVar("order_param_y", "The y component of the transformed order parameter");
  params.addCoupledVar("order_param_z", 0.0, "The z component of the transformed order parameter");
  params.addRequiredCoupledVar("f_q0", "The x component of the untransformed force");
  params.addRequiredCoupledVar("f_q1", "The y component of the untransformed force");
  params.addRequiredCoupledVar("f_q2", "The z component of the untransformed force");
  params.addRequiredCoupledVar("J_q0q0", "The xx component of the untransformed jacobian");
  params.addRequiredCoupledVar("J_q1q1", "The yy component of the untransformed jacobian");
  params.addRequiredCoupledVar("J_q2q2", "The zz component of the untransformed jacobian");
  params.addRequiredCoupledVar("J_q0q1", "The xy component of the untransformed jacobian");
  params.addRequiredCoupledVar("J_q1q2", "The yz component of the untransformed jacobian");
  params.addRequiredCoupledVar("J_q0q2", "The xz component of the untransformed jacobian");
  return params;
}

Transformed111KernelOp3::Transformed111KernelOp3(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _order_param_x_var(coupled("order_param_x")),
   _order_param_y_var(coupled("order_param_y")),
   _order_param_z_var(coupled("order_param_z")),
   _f_q0(coupledValue("f_q0")),
   _f_q1(coupledValue("f_q1")),
   _f_q2(coupledValue("f_q2")),
   _J_q0q0(coupledValue("J_q0q0")),
   _J_q1q1(coupledValue("J_q1q1")),
   _J_q2q2(coupledValue("J_q2q2")),
   _J_q0q1(coupledValue("J_q0q1")),
   _J_q1q2(coupledValue("J_q1q2")),
   _J_q0q2(coupledValue("J_q0q2"))
{
}

Real
Transformed111KernelOp3::computeQpResidual()

//
// TODO: Note that there is no reason this needs to be hardcoded, but this will be the first step. 
//       in general, this procedure should work for any transformation
//

//
//     calculate f' = S*f such that f'_(_component) is listed here.
//
{
  if (_component == 0)
  {
    return _test[_i][_qp] * (0.40824829046386301637*_f_q0[_qp] + 0.40824829046386301637*_f_q1[_qp] - 0.81649658092772603273*_f_q2[_qp]);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (-0.7071067811865475244*_f_q0[_qp] + 0.7071067811865475244*_f_q1[_qp]);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (0.57735026918962576451*_f_q0[_qp] + 0.57735026918962576451*_f_q1[_qp] + 0.57735026918962576451*_f_q2[_qp]);
  }
  else
    return 0.0;
}

Real
Transformed111KernelOp3::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667*(_J_q0q0[_qp] + 2.0*_J_q0q1[_qp] - 4.0*_J_q0q2[_qp] + _J_q1q1[_qp] - 4.0*_J_q1q2[_qp] + 4.0*_J_q2q2[_qp]));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.5*(_J_q0q0[_qp] - 2.0*_J_q0q1[_qp] + _J_q1q1[_qp]));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.33333333333333333333*(_J_q0q0[_qp] + 2.0*_J_q0q1[_qp] + 2.0*_J_q0q2[_qp] + _J_q1q1[_qp] + 2.0*_J_q1q2[_qp] + _J_q2q2[_qp]));
  }
  else
    return 0.0;
}


Real
Transformed111KernelOp3::computeQpOffDiagJacobian(unsigned int jvar)
{ 
  if (_component == 0)
  {
    if (jvar == _order_param_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.28867513459481288225*(-1.*_J_q0q0[_qp] + 2.*_J_q0q2[_qp] + _J_q1q1[_qp] - 2.*_J_q1q2[_qp]));
    }
    else if (jvar == _order_param_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.23570226039551584147*(_J_q0q0[_qp] + 2.*_J_q0q1[_qp] - 1.*_J_q0q2[_qp] + _J_q1q1[_qp] - 1.*_J_q1q2[_qp] - 2.*_J_q2q2[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _order_param_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.28867513459481288225*(-1.*_J_q0q0[_qp] + 2.*_J_q0q2[_qp] + _J_q1q1[_qp] - 2.*_J_q1q2[_qp]));
    }
    else if (jvar == _order_param_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.40824829046386301637*(-1.*_J_q0q0[_qp] - 1.*_J_q0q2[_qp] + _J_q1q1[_qp] + _J_q1q2[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _order_param_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.23570226039551584147*(_J_q0q0[_qp] + 2.*_J_q0q1[_qp] - 1.*_J_q0q2[_qp] + _J_q1q1[_qp] - 1.*_J_q1q2[_qp] - 2.*_J_q2q2[_qp]));
    }
    else if (jvar == _order_param_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.40824829046386301637*(-1.*_J_q0q0[_qp] - 1.*_J_q0q2[_qp] + _J_q1q1[_qp] + _J_q1q2[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
