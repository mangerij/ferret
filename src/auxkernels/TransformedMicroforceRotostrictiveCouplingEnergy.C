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

#include "TransformedMicroforceRotostrictiveCouplingEnergy.h"

registerMooseObject("FerretApp", TransformedMicroforceRotostrictiveCouplingEnergy);

InputParameters TransformedMicroforceRotostrictiveCouplingEnergy::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Computes the free energy density of the local electrostrictive coupling.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this auxkernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("u1_x", "The x component of the transformed displacement");
  params.addRequiredCoupledVar("u1_y", "The y component of the transformed displacement");
  params.addCoupledVar("u1_z", "The z component of the transformed displacement");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase vector");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase vector");
  params.addCoupledVar("antiphase_A_z", "The z component of the antiphase vector");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

TransformedMicroforceRotostrictiveCouplingEnergy::TransformedMicroforceRotostrictiveCouplingEnergy(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _u1_x_var(coupled("u1_x")),
   _u1_y_var(coupled("u1_y")),
   _u1_z_var(coupled("u1_z")),
   _antiphase_A_x_var(coupled("antiphase_A_x")),
   _antiphase_A_y_var(coupled("antiphase_A_y")),
   _antiphase_A_z_var(coupled("antiphase_A_z")),
   _u1_x_grad(coupledGradient("u1_x")),
   _u1_y_grad(coupledGradient("u1_y")),
   _u1_z_grad(coupledGradient("u1_z")),
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
   _r11(getMaterialProperty<Real>("r11")),
   _r12(getMaterialProperty<Real>("r12")),
   _r44(getMaterialProperty<Real>("r44")),
  _len_scale(getParam<Real>("len_scale"))
{
}

//qPPS^{-1}e'S simplifies... but minus sign and factor of 2 might need to be checked...

Real
TransformedMicroforceRotostrictiveCouplingEnergy::computeValue()
{
  if (_component == 0)
  {
   return 0.5*(-0.6666666666666666667*_antiphase_A_z[_qp]*_r44[_qp]*(4.*_u1_x_grad[_qp](0) - 3.4641016151377545871*(_u1_x_grad[_qp](1) + _u1_y_grad[_qp](0)) + 1.4142135623730950488*(_u1_x_grad[_qp](2) + _u1_z_grad[_qp](0)) +
      2.4494897427831780982*(_u1_y_grad[_qp](2) + _u1_z_grad[_qp](1)) - 4.*_u1_z_grad[_qp](2)) +
   1.3333333333333333333*_antiphase_A_y[_qp]*_r44[_qp]*(_u1_x_grad[_qp](0) - 3.*_u1_y_grad[_qp](1) + 1.4142135623730950488*(_u1_x_grad[_qp](2) + _u1_z_grad[_qp](0)) + 2.*_u1_z_grad[_qp](2)) +
   0.33333333333333333333*_antiphase_A_x[_qp]*_r11[_qp]*(_u1_x_grad[_qp](0) - 1.7320508075688772935*(_u1_x_grad[_qp](1) + _u1_y_grad[_qp](0)) + 3.*_u1_y_grad[_qp](1) +
      1.4142135623730950488*(_u1_x_grad[_qp](2) + _u1_z_grad[_qp](0)) - 2.4494897427831780982*(_u1_y_grad[_qp](2) + _u1_z_grad[_qp](1)) + 2.*_u1_z_grad[_qp](2)) +
   0.33333333333333333333*_antiphase_A_x[_qp]*_r12[_qp]*(5.*_u1_x_grad[_qp](0) + 1.7320508075688772935*(_u1_x_grad[_qp](1) + _u1_y_grad[_qp](0)) + 3.*_u1_y_grad[_qp](1) -
      1.4142135623730950488*(_u1_x_grad[_qp](2) + _u1_z_grad[_qp](0)) + 2.4494897427831780982*(_u1_y_grad[_qp](2) + _u1_z_grad[_qp](1)) + 4.*_u1_z_grad[_qp](2)));
  }
  else if (_component == 1)
  {
    return 0.5*(-0.6666666666666666667*_antiphase_A_z[_qp]*_r44[_qp]*(4.*_u1_x_grad[_qp](0) + 3.4641016151377545871*(_u1_x_grad[_qp](1) + _u1_y_grad[_qp](0)) + 1.4142135623730950488*(_u1_x_grad[_qp](2) + _u1_z_grad[_qp](0)) -
      2.4494897427831780982*(_u1_y_grad[_qp](2) + _u1_z_grad[_qp](1)) - 4.*_u1_z_grad[_qp](2)) +
   1.3333333333333333333*_antiphase_A_x[_qp]*_r44[_qp]*(_u1_x_grad[_qp](0) - 3.*_u1_y_grad[_qp](1) + 1.4142135623730950488*(_u1_x_grad[_qp](2) + _u1_z_grad[_qp](0)) + 2.*_u1_z_grad[_qp](2)) +
   0.33333333333333333333*_antiphase_A_y[_qp]*(_r11[_qp]*(_u1_x_grad[_qp](0) + 1.7320508075688772935*(_u1_x_grad[_qp](1) + _u1_y_grad[_qp](0)) + 3.*_u1_y_grad[_qp](1) +
         1.4142135623730950488*(_u1_x_grad[_qp](2) + _u1_z_grad[_qp](0)) + 2.4494897427831780982*(_u1_y_grad[_qp](2) + _u1_z_grad[_qp](1)) + 2.*_u1_z_grad[_qp](2)) +
      _r12[_qp]*(5.*_u1_x_grad[_qp](0) - 1.7320508075688772935*(_u1_x_grad[_qp](1) + _u1_y_grad[_qp](0)) + 3.*_u1_y_grad[_qp](1) - 1.4142135623730950488*(_u1_x_grad[_qp](2) + _u1_z_grad[_qp](0)) -
         2.4494897427831780982*(_u1_y_grad[_qp](2) + _u1_z_grad[_qp](1)) + 4.*_u1_z_grad[_qp](2))));
  }
  else if (_component == 2)
  {
    return 0.0;
  }
  else
    return 0.5*(-0.6666666666666666667*_antiphase_A_y[_qp]*_r44[_qp]*(4.*_u1_x_grad[_qp](0) + 3.4641016151377545871*(_u1_x_grad[_qp](1) + _u1_y_grad[_qp](0)) + 1.4142135623730950488*(_u1_x_grad[_qp](2) + _u1_z_grad[_qp](0)) -
      2.4494897427831780982*(_u1_y_grad[_qp](2) + _u1_z_grad[_qp](1)) - 4.*_u1_z_grad[_qp](2)) +
   0.6666666666666666667*_antiphase_A_x[_qp]*_r44[_qp]*(-4.*_u1_x_grad[_qp](0) + 3.4641016151377545871*(_u1_x_grad[_qp](1) + _u1_y_grad[_qp](0)) - 1.4142135623730950488*(_u1_x_grad[_qp](2) + _u1_z_grad[_qp](0)) -
      2.4494897427831780982*(_u1_y_grad[_qp](2) + _u1_z_grad[_qp](1)) + 4.*_u1_z_grad[_qp](2)) +
   0.6666666666666666667*_antiphase_A_z[_qp]*(_r11[_qp]*(2.*_u1_x_grad[_qp](0) - 1.4142135623730950488*(_u1_x_grad[_qp](2) + _u1_z_grad[_qp](0)) + _u1_z_grad[_qp](2)) +
      _r12[_qp]*(_u1_x_grad[_qp](0) + 3.*_u1_y_grad[_qp](1) + 1.4142135623730950488*(_u1_x_grad[_qp](2) + _u1_z_grad[_qp](0)) + 2.*_u1_z_grad[_qp](2))));
}
