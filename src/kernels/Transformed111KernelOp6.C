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

#include "Transformed111KernelOp6.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", Transformed111KernelOp6);

template<>
InputParameters validParams<Transformed111KernelOp6>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates the transformed residual for the local free energy which is an eighth order expansion in the polarization.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in. (0 for x, 1 for y, 2.0 for z)");
  params.addRequiredCoupledVar("order_param_x", "The x component of the transformed order parameter");
  params.addRequiredCoupledVar("order_param_y", "The y component of the transformed order parameter");
  params.addCoupledVar("order_param_z", 0.0, "The z component of the transformed order parameter");
  params.addRequiredCoupledVar("order_param2_x", "The x component of the second transformed order parameter");
  params.addRequiredCoupledVar("order_param2_y", "The y component of the second transformed order parameter");
  params.addCoupledVar("order_param2_z", 0.0, "The z component of the second transformed order parameter");
  params.addRequiredCoupledVar("f_q0", "The x component of the untransformed force");
  params.addRequiredCoupledVar("f_q1", "The y component of the untransformed force");
  params.addRequiredCoupledVar("f_q2", "The z component of the untransformed force");
  params.addRequiredCoupledVar("f_q3", "The x component of the untransformed force");
  params.addRequiredCoupledVar("f_q4", "The y component of the untransformed force");
  params.addRequiredCoupledVar("f_q5", "The z component of the untransformed force");
  params.addRequiredCoupledVar("J_q0q0", "");
  params.addRequiredCoupledVar("J_q1q1", "");
  params.addRequiredCoupledVar("J_q2q2", "");
  params.addRequiredCoupledVar("J_q3q3", "");
  params.addRequiredCoupledVar("J_q4q4", "");
  params.addRequiredCoupledVar("J_q5q5", "");
  params.addRequiredCoupledVar("J_q0q1", "");
  params.addRequiredCoupledVar("J_q1q2", "");
  params.addRequiredCoupledVar("J_q0q2", "");
  params.addRequiredCoupledVar("J_q0q3", "");
  params.addRequiredCoupledVar("J_q0q4", "");
  params.addRequiredCoupledVar("J_q0q5", "");
  params.addRequiredCoupledVar("J_q1q3", "");
  params.addRequiredCoupledVar("J_q1q4", "");
  params.addRequiredCoupledVar("J_q1q5", "");
  params.addRequiredCoupledVar("J_q2q3", "");
  params.addRequiredCoupledVar("J_q2q4", "");
  params.addRequiredCoupledVar("J_q2q5", "");
  params.addRequiredCoupledVar("J_q3q4", "");
  params.addRequiredCoupledVar("J_q3q5", "");
  params.addRequiredCoupledVar("J_q4q5", "");

  return params;
}

Transformed111KernelOp6::Transformed111KernelOp6(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _order_param_x_var(coupled("order_param_x")),
   _order_param_y_var(coupled("order_param_y")),
   _order_param_z_var(coupled("order_param_z")),
   _order_param2_x_var(coupled("order_param2_x")),
   _order_param2_y_var(coupled("order_param2_y")),
   _order_param2_z_var(coupled("order_param2_z")),
   _f_q0(coupledValue("f_q0")),
   _f_q1(coupledValue("f_q1")),
   _f_q2(coupledValue("f_q2")),
   _f_q3(coupledValue("f_q3")),
   _f_q4(coupledValue("f_q4")),
   _f_q5(coupledValue("f_q5")),
   _J_q0q0(coupledValue("J_q0q0")),
   _J_q1q1(coupledValue("J_q1q1")),
   _J_q2q2(coupledValue("J_q2q2")),
   _J_q3q3(coupledValue("J_q3q3")),
   _J_q4q4(coupledValue("J_q4q4")),
   _J_q5q5(coupledValue("J_q5q5")),
   _J_q0q1(coupledValue("J_q0q1")),
   _J_q1q2(coupledValue("J_q1q2")),
   _J_q0q2(coupledValue("J_q0q2")),
   _J_q0q3(coupledValue("J_q0q3")),
   _J_q0q4(coupledValue("J_q0q4")),
   _J_q0q5(coupledValue("J_q0q5")),
   _J_q1q3(coupledValue("J_q1q3")),
   _J_q1q4(coupledValue("J_q1q4")),
   _J_q1q5(coupledValue("J_q1q5")),
   _J_q2q3(coupledValue("J_q2q3")),
   _J_q2q4(coupledValue("J_q2q4")),
   _J_q2q5(coupledValue("J_q2q5")),
   _J_q3q4(coupledValue("J_q3q4")),
   _J_q3q5(coupledValue("J_q3q5")),
   _J_q4q5(coupledValue("J_q4q5"))
{
}

Real
Transformed111KernelOp6::computeQpResidual()
//
// TODO: Note that there is no reason this needs to be hardcoded, but this will be the first step. 
//       in general, this procedure should work for any transformation
//

//  index ordering like follows: f = f_Px,f_Py,f_Pz,f_Ax,f_Ay,f_Az
//
//  Note: switching o and o2 is possible due to the symmetricity of the coupling term.
//
//     calculate f' = S*f such that f'(_component) is listed here.

{
  if (_component == 0)
  {
    return _test[_i][_qp] * (0.40824829046386301637*_f_q0[_qp] - 0.7071067811865475244*_f_q1[_qp] + 0.57735026918962576451*_f_q2[_qp]);
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * (0.40824829046386301637*_f_q0[_qp] + 0.7071067811865475244*_f_q1[_qp] + 0.57735026918962576451*_f_q2[_qp]);
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * (-0.81649658092772603273*_f_q0[_qp] + 0.57735026918962576451*_f_q2[_qp]);
  }
  if (_component == 3)
  {
    return _test[_i][_qp] * (0.40824829046386301637*_f_q3[_qp] - 0.7071067811865475244*_f_q4[_qp] + 0.57735026918962576451*_f_q5[_qp]);
  }
  else if (_component == 4)
  {
    return _test[_i][_qp] * (0.40824829046386301637*_f_q3[_qp] + 0.7071067811865475244*_f_q4[_qp] + 0.57735026918962576451*_f_q5[_qp]);
  }
  else if (_component == 5)
  {
    return _test[_i][_qp] * (-0.81649658092772603273*_f_q3[_qp] + 0.57735026918962576451*_f_q5[_qp]);
  }
  else
    return 0.0;
}

Real
Transformed111KernelOp6::computeQpJacobian()
{
  if (_component == 0)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667*(_J_q0q0[_qp] - 3.4641016151377545871*_J_q0q1[_qp] + 
     2.8284271247461900976*_J_q0q2[_qp] + 3.0*_J_q1q1[_qp] - 4.8989794855663561964*_J_q1q2[_qp] + 2.0*_J_q2q2[_qp]));
  }
  else if (_component == 1)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667*(_J_q0q0[_qp] + 3.4641016151377545871*_J_q0q1[_qp] + 
     2.8284271247461900976*_J_q0q2[_qp] + 3.0*_J_q1q1[_qp] + 4.8989794855663561964*_J_q1q2[_qp] + 2.0*_J_q2q2[_qp]));
  }
  else if (_component == 2)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.33333333333333333333*(2.0*_J_q0q0[_qp] - 2.8284271247461900976*_J_q0q2[_qp] + _J_q2q2[_qp]));
  }
  else if (_component == 3)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667*(_J_q3q3[_qp] - 3.4641016151377545871*_J_q3q4[_qp] + 
     2.8284271247461900976*_J_q3q5[_qp] + 3.0*_J_q4q4[_qp] - 4.8989794855663561964*_J_q4q5[_qp] + 2.*_J_q5q5[_qp]));
  }
  else if (_component == 4)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667*(_J_q3q3[_qp] + 3.4641016151377545871*_J_q3q4[_qp] + 
     2.8284271247461900976*_J_q3q5[_qp] + 3.0*_J_q4q4[_qp] + 4.8989794855663561964*_J_q4q5[_qp] + 2.*_J_q5q5[_qp]));
  }
  else if (_component == 5)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (0.33333333333333333333*(2.0*_J_q3q3[_qp] - 2.8284271247461900976*_J_q3q5[_qp] + _J_q5q5[_qp]));
  }
  else 
    return 0.0;
}


Real
Transformed111KernelOp6::computeQpOffDiagJacobian(unsigned int jvar)
{ 
  if (_component == 0)
  {
    if (jvar == _order_param_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667*(_J_q0q0[_qp] + 2.8284271247461900976*_J_q0q2[_qp] - 3.*_J_q1q1[_qp] + 
     2.*_J_q2q2[_qp]));
    }
    else if (jvar == _order_param_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667*(-2.*_J_q0q0[_qp] + 3.4641016151377545871*_J_q0q1[_qp] - 
     1.4142135623730950488*_J_q0q2[_qp] - 2.4494897427831780982*_J_q1q2[_qp] + 2.*_J_q2q2[_qp]));
    }
    else if (jvar == _order_param2_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667*(_J_q0q3[_qp] - 1.7320508075688772935*_J_q0q4[_qp] + 
     1.4142135623730950488*_J_q0q5[_qp] - 1.7320508075688772935*_J_q1q3[_qp] + 3.*_J_q1q4[_qp] - 2.4494897427831780982*_J_q1q5[_qp] + 1.4142135623730950488*_J_q2q3[_qp] - 2.4494897427831780982*_J_q2q4[_qp] + 2.*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param2_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667*(_J_q0q3[_qp] + 1.7320508075688772935*_J_q0q4[_qp] + 
     1.4142135623730950488*_J_q0q5[_qp] - 1.7320508075688772935*_J_q1q3[_qp] - 3.*_J_q1q4[_qp] - 2.4494897427831780982*_J_q1q5[_qp] + 
     1.4142135623730950488*_J_q2q3[_qp] + 2.4494897427831780982*_J_q2q4[_qp] + 2.*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param2_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667*(-2.*_J_q0q3[_qp] + 1.4142135623730950488*_J_q0q5[_qp] + 
     3.4641016151377545871*_J_q1q3[_qp] - 2.4494897427831780982*_J_q1q5[_qp] - 2.8284271247461900976*_J_q2q3[_qp] + 2.0*_J_q2q5[_qp]));
    }
    else
      return 0.0;
  }
  else if (_component == 1)
  {
    if (jvar == _order_param_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667 * (_J_q0q0[_qp] + 2.8284271247461900976*_J_q0q2[_qp] - 3.0*_J_q1q1[_qp] + 
     2.*_J_q2q2[_qp]));
    }
    else if (jvar == _order_param_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667 * (-2.0*_J_q0q0[_qp] - 3.4641016151377545871*_J_q0q1[_qp] - 
     1.4142135623730950488*_J_q0q2[_qp] + 2.4494897427831780982*_J_q1q2[_qp] + 2.0*_J_q2q2[_qp]));
    }
    else if (jvar == _order_param2_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667 * (_J_q0q3[_qp] - 1.7320508075688772935*_J_q0q4[_qp] + 
     1.4142135623730950488*_J_q0q5[_qp] + 1.7320508075688772935*_J_q1q3[_qp] - 3.0*_J_q1q4[_qp] + 2.4494897427831780982*_J_q1q5[_qp] + 1.4142135623730950488*_J_q2q3[_qp] - 2.4494897427831780982*_J_q2q4[_qp] + 2.*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param2_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667 * (_J_q0q3[_qp] + 1.7320508075688772935*_J_q0q4[_qp] + 
     1.4142135623730950488*_J_q0q5[_qp] + 1.7320508075688772935*_J_q1q3[_qp] + 3.0*_J_q1q4[_qp] + 2.4494897427831780982*_J_q1q5[_qp] + 
     1.4142135623730950488*_J_q2q3[_qp] + 2.4494897427831780982*_J_q2q4[_qp] + 2.0*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param2_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667*(-2.0*_J_q0q3[_qp] + 1.4142135623730950488*_J_q0q5[_qp] - 
     3.4641016151377545871*_J_q1q3[_qp] + 2.4494897427831780982*_J_q1q5[_qp] - 2.8284271247461900976*_J_q2q3[_qp] + 2.*_J_q2q5[_qp]));
    }
    else
      return 0.0;
  }
  else if (_component == 2)
  {
    if (jvar == _order_param_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667 * (-2.0*_J_q0q0[_qp] + 3.4641016151377545871*_J_q0q1[_qp] - 
     1.4142135623730950488*_J_q0q2[_qp] - 2.4494897427831780982*_J_q1q2[_qp] + 2.0*_J_q2q2[_qp]));
    }
    else if (jvar == _order_param_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667 * (-2.0*_J_q0q0[_qp] - 3.4641016151377545871*_J_q0q1[_qp] - 
     1.4142135623730950488*_J_q0q2[_qp] + 2.4494897427831780982*_J_q1q2[_qp] + 2.0*_J_q2q2[_qp]));
    }
    else if (jvar == _order_param2_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667 * (-2.0*_J_q0q3[_qp] + 3.4641016151377545871*_J_q0q4[_qp] - 
     2.8284271247461900976*_J_q0q5[_qp] + 1.4142135623730950488*_J_q2q3[_qp] - 2.4494897427831780982*_J_q2q4[_qp] + 2.0*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param2_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.16666666666666666667*(-2.0*_J_q0q3[_qp] - 3.4641016151377545871*_J_q0q4[_qp] - 
     2.8284271247461900976*_J_q0q5[_qp] + 1.4142135623730950488*_J_q2q3[_qp] + 2.4494897427831780982*_J_q2q4[_qp] + 2.0*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param2_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] * (0.33333333333333333333*(2.0*_J_q0q3[_qp] - 1.4142135623730950488*_J_q0q5[_qp] - 
     1.4142135623730950488*_J_q2q3[_qp] + _J_q2q5[_qp]));
    }
    else
      return 0.0;
  }
  else if (_component == 3)
  {
    if (jvar == _order_param_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667* (_J_q0q3[_qp] - 1.7320508075688772935*_J_q0q4[_qp] + 
     1.4142135623730950488*_J_q0q5[_qp] - 1.7320508075688772935*_J_q1q3[_qp] + 3.*_J_q1q4[_qp] - 2.4494897427831780982*_J_q1q5[_qp] + 
     1.4142135623730950488*_J_q2q3[_qp] - 2.4494897427831780982*_J_q2q4[_qp] + 2.*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667* (_J_q0q3[_qp] - 1.7320508075688772935*_J_q0q4[_qp] + 
     1.4142135623730950488*_J_q0q5[_qp] + 1.7320508075688772935*_J_q1q3[_qp] - 3.*_J_q1q4[_qp] + 2.4494897427831780982*_J_q1q5[_qp] + 
     1.4142135623730950488*_J_q2q3[_qp] - 2.4494897427831780982*_J_q2q4[_qp] + 2.*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667*(-2.*_J_q0q3[_qp] + 3.4641016151377545871*_J_q0q4[_qp] - 
     2.8284271247461900976*_J_q0q5[_qp] + 1.4142135623730950488*_J_q2q3[_qp] - 2.4494897427831780982*_J_q2q4[_qp] + 2.*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param2_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667*(_J_q3q3[_qp] + 2.8284271247461900976*_J_q3q5[_qp] - 3.*_J_q4q4[_qp] + 
     2.*_J_q5q5[_qp]));
    }
    else if (jvar == _order_param2_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667*(-2.*_J_q3q3[_qp] + 3.4641016151377545871*_J_q3q4[_qp] - 
     1.4142135623730950488*_J_q3q5[_qp] - 
     2.4494897427831780982*_J_q4q5[_qp] + 2.*_J_q5q5[_qp]));
    }
    else
      return 0.0;
  }
  else if (_component == 4)
  {
    if (jvar == _order_param_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667*(_J_q0q3[_qp] + 1.7320508075688772935*_J_q0q4[_qp] + 
     1.4142135623730950488*_J_q0q5[_qp] - 1.7320508075688772935*_J_q1q3[_qp] - 3.*_J_q1q4[_qp] - 2.4494897427831780982*_J_q1q5[_qp] + 
     1.4142135623730950488*_J_q2q3[_qp] + 2.4494897427831780982*_J_q2q4[_qp] + 2.*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667* (_J_q0q3[_qp] + 1.7320508075688772935*_J_q0q4[_qp] + 
     1.4142135623730950488*_J_q0q5[_qp] + 1.7320508075688772935*_J_q1q3[_qp] + 3.*_J_q1q4[_qp] + 2.4494897427831780982*_J_q1q5[_qp] + 
     1.4142135623730950488*_J_q2q3[_qp] + 2.4494897427831780982*_J_q2q4[_qp] + 2.*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667*(-2.*_J_q0q3[_qp] - 3.4641016151377545871*_J_q0q4[_qp] - 
     2.8284271247461900976*_J_q0q5[_qp] + 1.4142135623730950488*_J_q2q3[_qp] + 2.4494897427831780982*_J_q2q4[_qp] + 2.*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param2_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667 * (_J_q3q3[_qp] + 2.8284271247461900976*_J_q3q5[_qp] - 3.*_J_q4q4[_qp] + 
     2.*_J_q5q5[_qp]));
    }
    else if (jvar == _order_param2_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667* (-2.*_J_q3q3[_qp] - 3.4641016151377545871*_J_q3q4[_qp] - 
     1.4142135623730950488*_J_q3q5[_qp] + 2.4494897427831780982*_J_q4q5[_qp] + 2.*_J_q5q5[_qp]));
    }
    else
      return 0.0;
  }
  else if (_component == 5)
  {
    if (jvar == _order_param_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667* (-2.*_J_q0q3[_qp] + 1.4142135623730950488*_J_q0q5[_qp] + 
     3.4641016151377545871*_J_q1q3[_qp] - 2.4494897427831780982*_J_q1q5[_qp] - 2.8284271247461900976*_J_q2q3[_qp] + 2.*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667* (-2.*_J_q0q3[_qp] + 1.4142135623730950488*_J_q0q5[_qp] - 
     3.4641016151377545871*_J_q1q3[_qp] + 2.4494897427831780982*_J_q1q5[_qp] - 2.8284271247461900976*_J_q2q3[_qp] + 2.*_J_q2q5[_qp]));
    }
    else if (jvar == _order_param_z_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.33333333333333333333* (2.*_J_q0q3[_qp] - 1.4142135623730950488*_J_q0q5[_qp] - 
     1.4142135623730950488*_J_q2q3[_qp] + _J_q2q5[_qp]));
    }
    else if (jvar == _order_param2_x_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667* (-2.*_J_q3q3[_qp] + 3.4641016151377545871*_J_q3q4[_qp] - 
     1.4142135623730950488*_J_q3q5[_qp] - 2.4494897427831780982*_J_q4q5[_qp] + 2.*_J_q5q5[_qp]));
    }
    else if (jvar == _order_param2_y_var)
    {
      return _test[_i][_qp] * _phi[_j][_qp] *  (0.16666666666666666667* (-2.*_J_q3q3[_qp] - 3.4641016151377545871*_J_q3q4[_qp] - 
     1.4142135623730950488*_J_q3q5[_qp] + 2.4494897427831780982*_J_q4q5[_qp] + 2.*_J_q5q5[_qp]));
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}
