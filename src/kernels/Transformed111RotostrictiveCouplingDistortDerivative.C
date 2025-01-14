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

#include "Transformed111RotostrictiveCouplingDistortDerivative.h"

class Transformed111RotostrictiveCouplingDistortDerivative;

registerMooseObject("FerretApp", Transformed111RotostrictiveCouplingDistortDerivative);

InputParameters Transformed111RotostrictiveCouplingDistortDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to the variation w.r.t antiphase of the rotostrictive coupling energy. Note: for cubic parent phase only only.");
  params.addRequiredCoupledVar("u_x", "The x component of the local elastic displacement");
  params.addRequiredCoupledVar("u_y", "The y component of the local elastic displacement");
  params.addCoupledVar("u_z", 0.0, "The z component of the local elastic displacement");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase tilt vector");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase tilt vector");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the antiphase tilt vector");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  return params;
}

Transformed111RotostrictiveCouplingDistortDerivative::Transformed111RotostrictiveCouplingDistortDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _u_x_var(coupled("u_x")),
   _u_y_var(coupled("u_y")),
   _u_z_var(coupled("u_z")),
   _antiphase_A_x_var(coupled("antiphase_A_x")),
   _antiphase_A_y_var(coupled("antiphase_A_y")),
   _antiphase_A_z_var(coupled("antiphase_A_z")),
   _u_x_grad(coupledGradient("u_x")),
   _u_y_grad(coupledGradient("u_y")),
   _u_z_grad(coupledGradient("u_z")),
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
   _r11(getMaterialProperty<Real>("r11")),
   _r12(getMaterialProperty<Real>("r12")),
   _r44(getMaterialProperty<Real>("r44"))
{
}

Real
Transformed111RotostrictiveCouplingDistortDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return -0.5*_test[_i][_qp] * (-0.33333333333333333333*(_antiphase_A_z[_qp]*(_r11[_qp]*(1.4142135623730950488*_u_x_grad[_qp](0) - 2.*_u_x_grad[_qp](2) - 1.4142135623730950488*_u_y_grad[_qp](1) - 2.*_u_z_grad[_qp](0)) -
         2.*_r44[_qp]*(1.4142135623730950488*_u_x_grad[_qp](0) + _u_x_grad[_qp](2) - 1.4142135623730950488*_u_y_grad[_qp](1) + _u_z_grad[_qp](0)) +
         _r12[_qp]*(-1.4142135623730950488*_u_x_grad[_qp](0) + 2.*_u_x_grad[_qp](2) + 1.4142135623730950488*_u_y_grad[_qp](1) + 2.*_u_z_grad[_qp](0))) +
      _antiphase_A_y[_qp]*(-1.*_r11[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) + _r12[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) +
         _r44[_qp]*(-4.*_u_x_grad[_qp](1) - 4.*_u_y_grad[_qp](0) + 2.8284271247461900976*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))))) +
   0.33333333333333333333*_antiphase_A_x[_qp]*(_r44[_qp]*(6.*_u_x_grad[_qp](0) + 2.8284271247461900976*_u_x_grad[_qp](2) - 2.*_u_y_grad[_qp](1) + 2.8284271247461900976*_u_z_grad[_qp](0) - 4.*_u_z_grad[_qp](2)) +
      _r11[_qp]*(3.*_u_x_grad[_qp](0) - 1.4142135623730950488*_u_x_grad[_qp](2) + _u_y_grad[_qp](1) - 1.4142135623730950488*_u_z_grad[_qp](0) + 2.*_u_z_grad[_qp](2)) +
      _r12[_qp]*(3.*_u_x_grad[_qp](0) + 1.4142135623730950488*_u_x_grad[_qp](2) + 5.*_u_y_grad[_qp](1) + 1.4142135623730950488*_u_z_grad[_qp](0) + 4.*_u_z_grad[_qp](2))));
  }
  else if (_component == 1)
  {
    return -0.5*_test[_i][_qp] * (0.33333333333333333333*_antiphase_A_z[_qp]*(1.4142135623730950488*_r11[_qp]*_u_x_grad[_qp](1) - 1.4142135623730950488*_r12[_qp]*_u_x_grad[_qp](1) + 1.4142135623730950488*_r11[_qp]*_u_y_grad[_qp](0) -
      1.4142135623730950488*_r12[_qp]*_u_y_grad[_qp](0) + 2.*_r11[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) - 2.*_r12[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) +
      2.*_r44[_qp]*(-1.4142135623730950488*_u_x_grad[_qp](1) - 1.4142135623730950488*_u_y_grad[_qp](0) + _u_y_grad[_qp](2) + _u_z_grad[_qp](1))) -
   0.33333333333333333333*_antiphase_A_x[_qp]*(-1.*_r11[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) +
      _r12[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) + _r44[_qp]*(-4.*_u_x_grad[_qp](1) - 4.*_u_y_grad[_qp](0) + 2.8284271247461900976*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)))) +
   0.33333333333333333333*_antiphase_A_y[_qp]*(-2.*_r44[_qp]*(_u_x_grad[_qp](0) + 1.4142135623730950488*_u_x_grad[_qp](2) - 3.*_u_y_grad[_qp](1) + 1.4142135623730950488*_u_z_grad[_qp](0) + 2.*_u_z_grad[_qp](2)) +
      _r11[_qp]*(_u_x_grad[_qp](0) + 1.4142135623730950488*_u_x_grad[_qp](2) + 3.*_u_y_grad[_qp](1) + 1.4142135623730950488*_u_z_grad[_qp](0) + 2.*_u_z_grad[_qp](2)) +
      _r12[_qp]*(5.*_u_x_grad[_qp](0) - 1.4142135623730950488*_u_x_grad[_qp](2) + 3.*_u_y_grad[_qp](1) - 1.4142135623730950488*_u_z_grad[_qp](0) + 4.*_u_z_grad[_qp](2))));
  }
  else if (_component == 2)
  {
    return -0.5*_test[_i][_qp] * (-0.33333333333333333333*_antiphase_A_x[_qp]*(_r11[_qp]*(1.4142135623730950488*_u_x_grad[_qp](0) - 2.*_u_x_grad[_qp](2) - 1.4142135623730950488*_u_y_grad[_qp](1) - 2.*_u_z_grad[_qp](0)) -
      2.*_r44[_qp]*(1.4142135623730950488*_u_x_grad[_qp](0) + _u_x_grad[_qp](2) - 1.4142135623730950488*_u_y_grad[_qp](1) + _u_z_grad[_qp](0)) +
      _r12[_qp]*(-1.4142135623730950488*_u_x_grad[_qp](0) + 2.*_u_x_grad[_qp](2) + 1.4142135623730950488*_u_y_grad[_qp](1) + 2.*_u_z_grad[_qp](0))) +
   0.33333333333333333333*_antiphase_A_y[_qp]*(1.4142135623730950488*_r11[_qp]*_u_x_grad[_qp](1) - 1.4142135623730950488*_r12[_qp]*_u_x_grad[_qp](1) + 1.4142135623730950488*_r11[_qp]*_u_y_grad[_qp](0) -
      1.4142135623730950488*_r12[_qp]*_u_y_grad[_qp](0) + 2.*_r11[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) - 2.*_r12[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) +
      2.*_r44[_qp]*(-1.4142135623730950488*_u_x_grad[_qp](1) - 1.4142135623730950488*_u_y_grad[_qp](0) + _u_y_grad[_qp](2) + _u_z_grad[_qp](1))) +
   0.6666666666666666667*_antiphase_A_z[_qp]*(-2.*_r44[_qp]*(_u_x_grad[_qp](0) + _u_y_grad[_qp](1) - 2.*_u_z_grad[_qp](2)) + _r11[_qp]*(_u_x_grad[_qp](0) + _u_y_grad[_qp](1) + _u_z_grad[_qp](2)) + 2.*_r12[_qp]*(_u_x_grad[_qp](0) + _u_y_grad[_qp](1) + _u_z_grad[_qp](2))));
  }
  else
    return 0.0;
}

Real
Transformed111RotostrictiveCouplingDistortDerivative::computeQpJacobian()
{
  if (_component == 0)
  {
    return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (0.33333333333333333333*(_r44[_qp]*(6.*_u_x_grad[_qp](0) + 2.8284271247461900976*_u_x_grad[_qp](2) - 2.*_u_y_grad[_qp](1) + 2.8284271247461900976*_u_z_grad[_qp](0) - 4.*_u_z_grad[_qp](2)) +
     _r11[_qp]*(3.*_u_x_grad[_qp](0) - 1.4142135623730950488*_u_x_grad[_qp](2) + _u_y_grad[_qp](1) - 1.4142135623730950488*_u_z_grad[_qp](0) + 2.*_u_z_grad[_qp](2)) +
     _r12[_qp]*(3.*_u_x_grad[_qp](0) + 1.4142135623730950488*_u_x_grad[_qp](2) + 5.*_u_y_grad[_qp](1) + 1.4142135623730950488*_u_z_grad[_qp](0) + 4.*_u_z_grad[_qp](2))));
  }
  else if (_component == 1)
  {
    return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (0.33333333333333333333*(-2.*_r44[_qp]*(_u_x_grad[_qp](0) + 1.4142135623730950488*_u_x_grad[_qp](2) - 3.*_u_y_grad[_qp](1) + 1.4142135623730950488*_u_z_grad[_qp](0) + 2.*_u_z_grad[_qp](2)) +
     _r11[_qp]*(_u_x_grad[_qp](0) + 1.4142135623730950488*_u_x_grad[_qp](2) + 3.*_u_y_grad[_qp](1) + 1.4142135623730950488*_u_z_grad[_qp](0) + 2.*_u_z_grad[_qp](2)) +
     _r12[_qp]*(5.*_u_x_grad[_qp](0) - 1.4142135623730950488*_u_x_grad[_qp](2) + 3.*_u_y_grad[_qp](1) - 1.4142135623730950488*_u_z_grad[_qp](0) + 4.*_u_z_grad[_qp](2))));
  }
  else if (_component == 2)
  {
    return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (-1.3333333333333333333*_r44[_qp]*(_u_x_grad[_qp](0) + _u_y_grad[_qp](1) - 2.*_u_z_grad[_qp](2)) + 0.66666666666666666667*_r11[_qp]*(_u_x_grad[_qp](0) + _u_y_grad[_qp](1) + _u_z_grad[_qp](2)) +
   1.3333333333333333333*_r12[_qp]*(_u_x_grad[_qp](0) + _u_y_grad[_qp](1) + _u_z_grad[_qp](2)));
  }
  else
    return 0.0;
}

Real
Transformed111RotostrictiveCouplingDistortDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiphase_A_y_var)
    {
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (_r44[_qp]*(1.3333333333333333333*_u_x_grad[_qp](1) + 1.3333333333333333333*_u_y_grad[_qp](0) - 0.9428090415820633659*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) +
   0.33333333333333333333*_r11[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) -
   0.33333333333333333333*_r12[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))));
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (_r12[_qp]*(0.4714045207910316829*_u_x_grad[_qp](0) - 0.6666666666666666667*_u_x_grad[_qp](2) - 0.4714045207910316829*_u_y_grad[_qp](1) - 0.6666666666666666667*_u_z_grad[_qp](0)) +
   _r11[_qp]*(-0.4714045207910316829*_u_x_grad[_qp](0) + 0.6666666666666666667*_u_x_grad[_qp](2) + 0.4714045207910316829*_u_y_grad[_qp](1) + 0.6666666666666666667*_u_z_grad[_qp](0)) +
   _r44[_qp]*(0.9428090415820633659*_u_x_grad[_qp](0) + 0.6666666666666666667*_u_x_grad[_qp](2) - 0.9428090415820633659*_u_y_grad[_qp](1) + 0.6666666666666666667*_u_z_grad[_qp](0)));
    }
    else if (jvar == _u_x_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](0)*_antiphase_A_z[_qp]*(-9.428090415820634*_r11[_qp] + 9.428090415820634*_r12[_qp] + 18.856180831641268*_r44[_qp]) +
   _grad_phi[_j][_qp](1)*_antiphase_A_y[_qp]*(6.666666666666666*_r11[_qp] - 6.666666666666666*_r12[_qp] + 26.666666666666664*_r44[_qp]) + _grad_phi[_j][_qp](0)*_antiphase_A_x[_qp]*(20.*_r11[_qp] + 20.*_r12[_qp] + 40.*_r44[_qp]) +
   _grad_phi[_j][_qp](2)*(_antiphase_A_z[_qp]*(13.333333333333332*_r11[_qp] - 13.333333333333332*_r12[_qp] + 13.333333333333332*_r44[_qp]) +
      _antiphase_A_x[_qp]*(-9.428090415820634*_r11[_qp] + 9.428090415820634*_r12[_qp] + 18.856180831641268*_r44[_qp])));
    }
    else if (jvar == _u_y_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](2)*_antiphase_A_y[_qp]*(9.428090415820634*_r11[_qp] - 9.428090415820634*_r12[_qp] - 18.856180831641268*_r44[_qp]) +
   _grad_phi[_j][_qp](1)*_antiphase_A_z[_qp]*(9.428090415820634*_r11[_qp] - 9.428090415820634*_r12[_qp] - 18.856180831641268*_r44[_qp]) +
   _grad_phi[_j][_qp](1)*_antiphase_A_x[_qp]*(6.666666666666667*_r11[_qp] + 33.333333333333336*_r12[_qp] - 13.333333333333334*_r44[_qp]) +
   _grad_phi[_j][_qp](0)*_antiphase_A_y[_qp]*(6.666666666666667*_r11[_qp] - 6.666666666666667*_r12[_qp] + 26.666666666666668*_r44[_qp]));
    }
    else if (jvar == _u_z_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](2)*_antiphase_A_x[_qp]*(13.333333333333332*_r11[_qp] + 26.666666666666664*_r12[_qp] - 26.666666666666664*_r44[_qp]) +
   _grad_phi[_j][_qp](1)*_antiphase_A_y[_qp]*(9.428090415820634*_r11[_qp] - 9.428090415820634*_r12[_qp] - 18.856180831641268*_r44[_qp]) +
   _grad_phi[_j][_qp](0)*_antiphase_A_z[_qp]*(13.333333333333332*_r11[_qp] - 13.333333333333332*_r12[_qp] + 13.333333333333332*_r44[_qp]) +
   _grad_phi[_j][_qp](0)*_antiphase_A_x[_qp]*(-9.428090415820634*_r11[_qp] + 9.428090415820634*_r12[_qp] + 18.856180831641268*_r44[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 1)
  {
    if (jvar == _antiphase_A_x_var)
    {
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (_r44[_qp]*(1.3333333333333333333*_u_x_grad[_qp](1) + 1.3333333333333333333*_u_y_grad[_qp](0) - 0.9428090415820633659*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) +
   0.33333333333333333333*_r11[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) -
   0.33333333333333333333*_r12[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))));
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (0.47140452079103168293*_r11[_qp]*_u_x_grad[_qp](1) - 0.4714045207910316829*_r12[_qp]*_u_x_grad[_qp](1) + 0.47140452079103168293*_r11[_qp]*_u_y_grad[_qp](0) - 0.4714045207910316829*_r12[_qp]*_u_y_grad[_qp](0) +
   0.6666666666666666667*_r11[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) - 0.6666666666666666667*_r12[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) +
   0.6666666666666666667*_r44[_qp]*(-1.4142135623730950488*_u_x_grad[_qp](1) - 1.4142135623730950488*_u_y_grad[_qp](0) + _u_y_grad[_qp](2) + _u_z_grad[_qp](1)));
    }
    else if (jvar == _u_x_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](2)*_antiphase_A_y[_qp]*(9.428090415820634*_r11[_qp] - 9.428090415820634*_r12[_qp] - 18.856180831641268*_r44[_qp]) +
   _grad_phi[_j][_qp](1)*_antiphase_A_z[_qp]*(9.428090415820634*_r11[_qp] - 9.428090415820634*_r12[_qp] - 18.856180831641268*_r44[_qp]) +
   _grad_phi[_j][_qp](0)*_antiphase_A_y[_qp]*(6.666666666666667*_r11[_qp] + 33.333333333333336*_r12[_qp] - 13.333333333333334*_r44[_qp]) +
   _grad_phi[_j][_qp](1)*_antiphase_A_x[_qp]*(6.666666666666667*_r11[_qp] - 6.666666666666667*_r12[_qp] + 26.666666666666668*_r44[_qp]));
    }
    else if (jvar == _u_y_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](0)*_antiphase_A_z[_qp]*(9.428090415820634*_r11[_qp] - 9.428090415820634*_r12[_qp] - 18.856180831641268*_r44[_qp]) +
   _grad_phi[_j][_qp](0)*_antiphase_A_x[_qp]*(6.666666666666666*_r11[_qp] - 6.666666666666666*_r12[_qp] + 26.666666666666664*_r44[_qp]) + _grad_phi[_j][_qp](1)*_antiphase_A_y[_qp]*(20.*_r11[_qp] + 20.*_r12[_qp] + 40.*_r44[_qp]) +
   _grad_phi[_j][_qp](2)*(_antiphase_A_x[_qp]*(9.428090415820634*_r11[_qp] - 9.428090415820634*_r12[_qp] - 18.856180831641268*_r44[_qp]) +
      _antiphase_A_z[_qp]*(13.333333333333332*_r11[_qp] - 13.333333333333332*_r12[_qp] + 13.333333333333332*_r44[_qp])));
    }
    else if (jvar == _u_z_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](2)*_antiphase_A_y[_qp]*(13.333333333333334*_r11[_qp] + 26.666666666666668*_r12[_qp] - 26.666666666666668*_r44[_qp]) +
   _grad_phi[_j][_qp](1)*_antiphase_A_x[_qp]*(9.428090415820634*_r11[_qp] - 9.428090415820634*_r12[_qp] - 18.856180831641268*_r44[_qp]) +
   _grad_phi[_j][_qp](0)*_antiphase_A_y[_qp]*(9.428090415820634*_r11[_qp] - 9.428090415820634*_r12[_qp] - 18.856180831641268*_r44[_qp]) +
   _grad_phi[_j][_qp](1)*_antiphase_A_z[_qp]*(13.333333333333334*_r11[_qp] - 13.333333333333334*_r12[_qp] + 13.333333333333334*_r44[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else if (_component == 2)
  {
    if (jvar == _antiphase_A_x_var)
    {
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (_r12[_qp]*(0.4714045207910316829*_u_x_grad[_qp](0) - 0.6666666666666666667*_u_x_grad[_qp](2) - 0.4714045207910316829*_u_y_grad[_qp](1) - 0.6666666666666666667*_u_z_grad[_qp](0)) +
   _r11[_qp]*(-0.4714045207910316829*_u_x_grad[_qp](0) + 0.6666666666666666667*_u_x_grad[_qp](2) + 0.4714045207910316829*_u_y_grad[_qp](1) + 0.6666666666666666667*_u_z_grad[_qp](0)) +
   _r44[_qp]*(0.9428090415820633659*_u_x_grad[_qp](0) + 0.6666666666666666667*_u_x_grad[_qp](2) - 0.9428090415820633659*_u_y_grad[_qp](1) + 0.6666666666666666667*_u_z_grad[_qp](0)));
    }
    else if (jvar == _antiphase_A_y_var)
    {
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (0.47140452079103168293*_r11[_qp]*_u_x_grad[_qp](1) - 0.4714045207910316829*_r12[_qp]*_u_x_grad[_qp](1) + 0.47140452079103168293*_r11[_qp]*_u_y_grad[_qp](0) - 0.4714045207910316829*_r12[_qp]*_u_y_grad[_qp](0) +
   0.6666666666666666667*_r11[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) - 0.6666666666666666667*_r12[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) +
   0.6666666666666666667*_r44[_qp]*(-1.4142135623730950488*_u_x_grad[_qp](1) - 1.4142135623730950488*_u_y_grad[_qp](0) + _u_y_grad[_qp](2) + _u_z_grad[_qp](1)));
    }
    else if (jvar == _u_x_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](0)*_antiphase_A_z[_qp]*(13.333333333333334*_r11[_qp] + 26.666666666666668*_r12[_qp] - 26.666666666666668*_r44[_qp]) +
   _grad_phi[_j][_qp](1)*_antiphase_A_y[_qp]*(9.428090415820634*_r11[_qp] - 9.428090415820634*_r12[_qp] - 18.856180831641268*_r44[_qp]) +
   _grad_phi[_j][_qp](2)*_antiphase_A_x[_qp]*(13.333333333333334*_r11[_qp] - 13.333333333333334*_r12[_qp] + 13.333333333333334*_r44[_qp]) +
   _grad_phi[_j][_qp](0)*_antiphase_A_x[_qp]*(-9.428090415820634*_r11[_qp] + 9.428090415820634*_r12[_qp] + 18.856180831641268*_r44[_qp]));
    }
    else if (jvar == _u_y_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](1)*_antiphase_A_z[_qp]*(13.333333333333334*_r11[_qp] + 26.666666666666668*_r12[_qp] - 26.666666666666668*_r44[_qp]) +
   _grad_phi[_j][_qp](1)*_antiphase_A_x[_qp]*(9.428090415820634*_r11[_qp] - 9.428090415820634*_r12[_qp] - 18.856180831641268*_r44[_qp]) +
   _grad_phi[_j][_qp](0)*_antiphase_A_y[_qp]*(9.428090415820634*_r11[_qp] - 9.428090415820634*_r12[_qp] - 18.856180831641268*_r44[_qp]) +
   _grad_phi[_j][_qp](2)*_antiphase_A_y[_qp]*(13.333333333333334*_r11[_qp] - 13.333333333333334*_r12[_qp] + 13.333333333333334*_r44[_qp]));
    }
    else if (jvar == _u_z_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](0)*_antiphase_A_x[_qp]*(13.333333333333334*_r11[_qp] - 13.333333333333334*_r12[_qp] + 13.333333333333334*_r44[_qp]) +
   _grad_phi[_j][_qp](1)*_antiphase_A_y[_qp]*(13.333333333333334*_r11[_qp] - 13.333333333333334*_r12[_qp] + 13.333333333333334*_r44[_qp]) +
   _grad_phi[_j][_qp](2)*_antiphase_A_z[_qp]*(13.333333333333334*_r11[_qp] + 26.666666666666668*_r12[_qp] + 53.333333333333336*_r44[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
