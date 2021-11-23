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

#include "Transformed111ElectrostrictiveCouplingPolarDerivative.h"

class Transformed111ElectrostrictiveCouplingPolarDerivative;

registerMooseObject("FerretApp", Transformed111ElectrostrictiveCouplingPolarDerivative);

InputParameters Transformed111ElectrostrictiveCouplingPolarDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to the variation w.r.t polarization of the electrostrictive coupling energy. Note: for cubic parent phase only only.");
  params.addRequiredCoupledVar("u_x", "The x component of the local elastic displacement");
  params.addRequiredCoupledVar("u_y", "The y component of the local elastic displacement");
  params.addCoupledVar("u_z", 0.0, "The z component of the local elastic displacement");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  return params;
}

Transformed111ElectrostrictiveCouplingPolarDerivative::Transformed111ElectrostrictiveCouplingPolarDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _u_x_var(coupled("u_x")),
   _u_y_var(coupled("u_y")),
   _u_z_var(coupled("u_z")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _u_x_grad(coupledGradient("u_x")),
   _u_y_grad(coupledGradient("u_y")),
   _u_z_grad(coupledGradient("u_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _q11(getMaterialProperty<Real>("q11")),
   _q12(getMaterialProperty<Real>("q12")),
   _q44(getMaterialProperty<Real>("q44"))
{
}

Real
Transformed111ElectrostrictiveCouplingPolarDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return -0.5*_test[_i][_qp] * (-0.33333333333333333333*(_polar_z[_qp]*(_q11[_qp]*(1.4142135623730950488*_u_x_grad[_qp](0) - 2.*_u_x_grad[_qp](2) - 1.4142135623730950488*_u_y_grad[_qp](1) - 2.*_u_z_grad[_qp](0)) - 
         2.*_q44[_qp]*(1.4142135623730950488*_u_x_grad[_qp](0) + _u_x_grad[_qp](2) - 1.4142135623730950488*_u_y_grad[_qp](1) + _u_z_grad[_qp](0)) + 
         _q12[_qp]*(-1.4142135623730950488*_u_x_grad[_qp](0) + 2.*_u_x_grad[_qp](2) + 1.4142135623730950488*_u_y_grad[_qp](1) + 2.*_u_z_grad[_qp](0))) + 
      _polar_y[_qp]*(-1.*_q11[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) + _q12[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) + 
         _q44[_qp]*(-4.*_u_x_grad[_qp](1) - 4.*_u_y_grad[_qp](0) + 2.8284271247461900976*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))))) + 
   0.33333333333333333333*_polar_x[_qp]*(_q44[_qp]*(6.*_u_x_grad[_qp](0) + 2.8284271247461900976*_u_x_grad[_qp](2) - 2.*_u_y_grad[_qp](1) + 2.8284271247461900976*_u_z_grad[_qp](0) - 4.*_u_z_grad[_qp](2)) + 
      _q11[_qp]*(3.*_u_x_grad[_qp](0) - 1.4142135623730950488*_u_x_grad[_qp](2) + _u_y_grad[_qp](1) - 1.4142135623730950488*_u_z_grad[_qp](0) + 2.*_u_z_grad[_qp](2)) + 
      _q12[_qp]*(3.*_u_x_grad[_qp](0) + 1.4142135623730950488*_u_x_grad[_qp](2) + 5.*_u_y_grad[_qp](1) + 1.4142135623730950488*_u_z_grad[_qp](0) + 4.*_u_z_grad[_qp](2))));
  }
  else if (_component == 1)
  {
    return -0.5*_test[_i][_qp] * (0.33333333333333333333*_polar_z[_qp]*(1.4142135623730950488*_q11[_qp]*_u_x_grad[_qp](1) - 1.4142135623730950488*_q12[_qp]*_u_x_grad[_qp](1) + 1.4142135623730950488*_q11[_qp]*_u_y_grad[_qp](0) - 
      1.4142135623730950488*_q12[_qp]*_u_y_grad[_qp](0) + 2.*_q11[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) - 2.*_q12[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) + 
      2.*_q44[_qp]*(-1.4142135623730950488*_u_x_grad[_qp](1) - 1.4142135623730950488*_u_y_grad[_qp](0) + _u_y_grad[_qp](2) + _u_z_grad[_qp](1))) - 
   0.33333333333333333333*_polar_x[_qp]*(-1.*_q11[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) + 
      _q12[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) + _q44[_qp]*(-4.*_u_x_grad[_qp](1) - 4.*_u_y_grad[_qp](0) + 2.8284271247461900976*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)))) + 
   0.33333333333333333333*_polar_y[_qp]*(-2.*_q44[_qp]*(_u_x_grad[_qp](0) + 1.4142135623730950488*_u_x_grad[_qp](2) - 3.*_u_y_grad[_qp](1) + 1.4142135623730950488*_u_z_grad[_qp](0) + 2.*_u_z_grad[_qp](2)) + 
      _q11[_qp]*(_u_x_grad[_qp](0) + 1.4142135623730950488*_u_x_grad[_qp](2) + 3.*_u_y_grad[_qp](1) + 1.4142135623730950488*_u_z_grad[_qp](0) + 2.*_u_z_grad[_qp](2)) + 
      _q12[_qp]*(5.*_u_x_grad[_qp](0) - 1.4142135623730950488*_u_x_grad[_qp](2) + 3.*_u_y_grad[_qp](1) - 1.4142135623730950488*_u_z_grad[_qp](0) + 4.*_u_z_grad[_qp](2))));
  }
  else if (_component == 2)
  {
    return -0.5*_test[_i][_qp] * (-0.33333333333333333333*_polar_x[_qp]*(_q11[_qp]*(1.4142135623730950488*_u_x_grad[_qp](0) - 2.*_u_x_grad[_qp](2) - 1.4142135623730950488*_u_y_grad[_qp](1) - 2.*_u_z_grad[_qp](0)) - 
      2.*_q44[_qp]*(1.4142135623730950488*_u_x_grad[_qp](0) + _u_x_grad[_qp](2) - 1.4142135623730950488*_u_y_grad[_qp](1) + _u_z_grad[_qp](0)) + 
      _q12[_qp]*(-1.4142135623730950488*_u_x_grad[_qp](0) + 2.*_u_x_grad[_qp](2) + 1.4142135623730950488*_u_y_grad[_qp](1) + 2.*_u_z_grad[_qp](0))) + 
   0.33333333333333333333*_polar_y[_qp]*(1.4142135623730950488*_q11[_qp]*_u_x_grad[_qp](1) - 1.4142135623730950488*_q12[_qp]*_u_x_grad[_qp](1) + 1.4142135623730950488*_q11[_qp]*_u_y_grad[_qp](0) - 
      1.4142135623730950488*_q12[_qp]*_u_y_grad[_qp](0) + 2.*_q11[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) - 2.*_q12[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) + 
      2.*_q44[_qp]*(-1.4142135623730950488*_u_x_grad[_qp](1) - 1.4142135623730950488*_u_y_grad[_qp](0) + _u_y_grad[_qp](2) + _u_z_grad[_qp](1))) + 
   0.6666666666666666667*_polar_z[_qp]*(-2.*_q44[_qp]*(_u_x_grad[_qp](0) + _u_y_grad[_qp](1) - 2.*_u_z_grad[_qp](2)) + _q11[_qp]*(_u_x_grad[_qp](0) + _u_y_grad[_qp](1) + _u_z_grad[_qp](2)) + 2.*_q12[_qp]*(_u_x_grad[_qp](0) + _u_y_grad[_qp](1) + _u_z_grad[_qp](2))));
  }
  else
    return 0.0;
}

Real
Transformed111ElectrostrictiveCouplingPolarDerivative::computeQpJacobian()
{
  if (_component == 0)
  {
    return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (0.33333333333333333333*(_q44[_qp]*(6.*_u_x_grad[_qp](0) + 2.8284271247461900976*_u_x_grad[_qp](2) - 2.*_u_y_grad[_qp](1) + 2.8284271247461900976*_u_z_grad[_qp](0) - 4.*_u_z_grad[_qp](2)) + 
     _q11[_qp]*(3.*_u_x_grad[_qp](0) - 1.4142135623730950488*_u_x_grad[_qp](2) + _u_y_grad[_qp](1) - 1.4142135623730950488*_u_z_grad[_qp](0) + 2.*_u_z_grad[_qp](2)) + 
     _q12[_qp]*(3.*_u_x_grad[_qp](0) + 1.4142135623730950488*_u_x_grad[_qp](2) + 5.*_u_y_grad[_qp](1) + 1.4142135623730950488*_u_z_grad[_qp](0) + 4.*_u_z_grad[_qp](2))));
  }
  else if (_component == 1)
  {
    return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (0.33333333333333333333*(-2.*_q44[_qp]*(_u_x_grad[_qp](0) + 1.4142135623730950488*_u_x_grad[_qp](2) - 3.*_u_y_grad[_qp](1) + 1.4142135623730950488*_u_z_grad[_qp](0) + 2.*_u_z_grad[_qp](2)) + 
     _q11[_qp]*(_u_x_grad[_qp](0) + 1.4142135623730950488*_u_x_grad[_qp](2) + 3.*_u_y_grad[_qp](1) + 1.4142135623730950488*_u_z_grad[_qp](0) + 2.*_u_z_grad[_qp](2)) + 
     _q12[_qp]*(5.*_u_x_grad[_qp](0) - 1.4142135623730950488*_u_x_grad[_qp](2) + 3.*_u_y_grad[_qp](1) - 1.4142135623730950488*_u_z_grad[_qp](0) + 4.*_u_z_grad[_qp](2))));
  }
  else if (_component == 2)
  {
    return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (-1.3333333333333333333*_q44[_qp]*(_u_x_grad[_qp](0) + _u_y_grad[_qp](1) - 2.*_u_z_grad[_qp](2)) + 0.66666666666666666667*_q11[_qp]*(_u_x_grad[_qp](0) + _u_y_grad[_qp](1) + _u_z_grad[_qp](2)) + 
   1.3333333333333333333*_q12[_qp]*(_u_x_grad[_qp](0) + _u_y_grad[_qp](1) + _u_z_grad[_qp](2)));
  }
  else
    return 0.0;
}

Real
Transformed111ElectrostrictiveCouplingPolarDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _polar_y_var)
    {
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (_q44[_qp]*(1.3333333333333333333*_u_x_grad[_qp](1) + 1.3333333333333333333*_u_y_grad[_qp](0) - 0.9428090415820633659*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) + 
   0.33333333333333333333*_q11[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) - 
   0.33333333333333333333*_q12[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))));
    }
    else if (jvar == _polar_z_var)
    {
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (_q12[_qp]*(0.4714045207910316829*_u_x_grad[_qp](0) - 0.6666666666666666667*_u_x_grad[_qp](2) - 0.4714045207910316829*_u_y_grad[_qp](1) - 0.6666666666666666667*_u_z_grad[_qp](0)) + 
   _q11[_qp]*(-0.4714045207910316829*_u_x_grad[_qp](0) + 0.6666666666666666667*_u_x_grad[_qp](2) + 0.4714045207910316829*_u_y_grad[_qp](1) + 0.6666666666666666667*_u_z_grad[_qp](0)) + 
   _q44[_qp]*(0.9428090415820633659*_u_x_grad[_qp](0) + 0.6666666666666666667*_u_x_grad[_qp](2) - 0.9428090415820633659*_u_y_grad[_qp](1) + 0.6666666666666666667*_u_z_grad[_qp](0)));
    }
    else if (jvar == _u_x_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](0)*_polar_z[_qp]*(-9.428090415820634*_q11[_qp] + 9.428090415820634*_q12[_qp] + 18.856180831641268*_q44[_qp]) + 
   _grad_phi[_j][_qp](1)*_polar_y[_qp]*(6.666666666666666*_q11[_qp] - 6.666666666666666*_q12[_qp] + 26.666666666666664*_q44[_qp]) + _grad_phi[_j][_qp](0)*_polar_x[_qp]*(20.*_q11[_qp] + 20.*_q12[_qp] + 40.*_q44[_qp]) + 
   _grad_phi[_j][_qp](2)*(_polar_z[_qp]*(13.333333333333332*_q11[_qp] - 13.333333333333332*_q12[_qp] + 13.333333333333332*_q44[_qp]) + 
      _polar_x[_qp]*(-9.428090415820634*_q11[_qp] + 9.428090415820634*_q12[_qp] + 18.856180831641268*_q44[_qp])));
    }
    else if (jvar == _u_y_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](2)*_polar_y[_qp]*(9.428090415820634*_q11[_qp] - 9.428090415820634*_q12[_qp] - 18.856180831641268*_q44[_qp]) + 
   _grad_phi[_j][_qp](1)*_polar_z[_qp]*(9.428090415820634*_q11[_qp] - 9.428090415820634*_q12[_qp] - 18.856180831641268*_q44[_qp]) + 
   _grad_phi[_j][_qp](1)*_polar_x[_qp]*(6.666666666666667*_q11[_qp] + 33.333333333333336*_q12[_qp] - 13.333333333333334*_q44[_qp]) + 
   _grad_phi[_j][_qp](0)*_polar_y[_qp]*(6.666666666666667*_q11[_qp] - 6.666666666666667*_q12[_qp] + 26.666666666666668*_q44[_qp]));
    }
    else if (jvar == _u_z_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](2)*_polar_x[_qp]*(13.333333333333332*_q11[_qp] + 26.666666666666664*_q12[_qp] - 26.666666666666664*_q44[_qp]) + 
   _grad_phi[_j][_qp](1)*_polar_y[_qp]*(9.428090415820634*_q11[_qp] - 9.428090415820634*_q12[_qp] - 18.856180831641268*_q44[_qp]) + 
   _grad_phi[_j][_qp](0)*_polar_z[_qp]*(13.333333333333332*_q11[_qp] - 13.333333333333332*_q12[_qp] + 13.333333333333332*_q44[_qp]) + 
   _grad_phi[_j][_qp](0)*_polar_x[_qp]*(-9.428090415820634*_q11[_qp] + 9.428090415820634*_q12[_qp] + 18.856180831641268*_q44[_qp]));
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
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (_q44[_qp]*(1.3333333333333333333*_u_x_grad[_qp](1) + 1.3333333333333333333*_u_y_grad[_qp](0) - 0.9428090415820633659*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) + 
   0.33333333333333333333*_q11[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))) - 
   0.33333333333333333333*_q12[_qp]*(_u_x_grad[_qp](1) + _u_y_grad[_qp](0) + 1.4142135623730950488*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1))));
    }
    else if (jvar == _polar_z_var)
    {
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (0.47140452079103168293*_q11[_qp]*_u_x_grad[_qp](1) - 0.4714045207910316829*_q12[_qp]*_u_x_grad[_qp](1) + 0.47140452079103168293*_q11[_qp]*_u_y_grad[_qp](0) - 0.4714045207910316829*_q12[_qp]*_u_y_grad[_qp](0) + 
   0.6666666666666666667*_q11[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) - 0.6666666666666666667*_q12[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) + 
   0.6666666666666666667*_q44[_qp]*(-1.4142135623730950488*_u_x_grad[_qp](1) - 1.4142135623730950488*_u_y_grad[_qp](0) + _u_y_grad[_qp](2) + _u_z_grad[_qp](1)));
    }
    else if (jvar == _u_x_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](2)*_polar_y[_qp]*(9.428090415820634*_q11[_qp] - 9.428090415820634*_q12[_qp] - 18.856180831641268*_q44[_qp]) + 
   _grad_phi[_j][_qp](1)*_polar_z[_qp]*(9.428090415820634*_q11[_qp] - 9.428090415820634*_q12[_qp] - 18.856180831641268*_q44[_qp]) + 
   _grad_phi[_j][_qp](0)*_polar_y[_qp]*(6.666666666666667*_q11[_qp] + 33.333333333333336*_q12[_qp] - 13.333333333333334*_q44[_qp]) + 
   _grad_phi[_j][_qp](1)*_polar_x[_qp]*(6.666666666666667*_q11[_qp] - 6.666666666666667*_q12[_qp] + 26.666666666666668*_q44[_qp]));
    }
    else if (jvar == _u_y_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](0)*_polar_z[_qp]*(9.428090415820634*_q11[_qp] - 9.428090415820634*_q12[_qp] - 18.856180831641268*_q44[_qp]) + 
   _grad_phi[_j][_qp](0)*_polar_x[_qp]*(6.666666666666666*_q11[_qp] - 6.666666666666666*_q12[_qp] + 26.666666666666664*_q44[_qp]) + _grad_phi[_j][_qp](1)*_polar_y[_qp]*(20.*_q11[_qp] + 20.*_q12[_qp] + 40.*_q44[_qp]) + 
   _grad_phi[_j][_qp](2)*(_polar_x[_qp]*(9.428090415820634*_q11[_qp] - 9.428090415820634*_q12[_qp] - 18.856180831641268*_q44[_qp]) + 
      _polar_z[_qp]*(13.333333333333332*_q11[_qp] - 13.333333333333332*_q12[_qp] + 13.333333333333332*_q44[_qp])));
    }
    else if (jvar == _u_z_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](2)*_polar_y[_qp]*(13.333333333333334*_q11[_qp] + 26.666666666666668*_q12[_qp] - 26.666666666666668*_q44[_qp]) + 
   _grad_phi[_j][_qp](1)*_polar_x[_qp]*(9.428090415820634*_q11[_qp] - 9.428090415820634*_q12[_qp] - 18.856180831641268*_q44[_qp]) + 
   _grad_phi[_j][_qp](0)*_polar_y[_qp]*(9.428090415820634*_q11[_qp] - 9.428090415820634*_q12[_qp] - 18.856180831641268*_q44[_qp]) + 
   _grad_phi[_j][_qp](1)*_polar_z[_qp]*(13.333333333333334*_q11[_qp] - 13.333333333333334*_q12[_qp] + 13.333333333333334*_q44[_qp]));
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
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (_q12[_qp]*(0.4714045207910316829*_u_x_grad[_qp](0) - 0.6666666666666666667*_u_x_grad[_qp](2) - 0.4714045207910316829*_u_y_grad[_qp](1) - 0.6666666666666666667*_u_z_grad[_qp](0)) + 
   _q11[_qp]*(-0.4714045207910316829*_u_x_grad[_qp](0) + 0.6666666666666666667*_u_x_grad[_qp](2) + 0.4714045207910316829*_u_y_grad[_qp](1) + 0.6666666666666666667*_u_z_grad[_qp](0)) + 
   _q44[_qp]*(0.9428090415820633659*_u_x_grad[_qp](0) + 0.6666666666666666667*_u_x_grad[_qp](2) - 0.9428090415820633659*_u_y_grad[_qp](1) + 0.6666666666666666667*_u_z_grad[_qp](0)));
    }
    else if (jvar == _polar_y_var)
    {
      return -0.5*_test[_i][_qp] * _phi[_j][_qp] * (0.47140452079103168293*_q11[_qp]*_u_x_grad[_qp](1) - 0.4714045207910316829*_q12[_qp]*_u_x_grad[_qp](1) + 0.47140452079103168293*_q11[_qp]*_u_y_grad[_qp](0) - 0.4714045207910316829*_q12[_qp]*_u_y_grad[_qp](0) + 
   0.6666666666666666667*_q11[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) - 0.6666666666666666667*_q12[_qp]*(_u_y_grad[_qp](2) + _u_z_grad[_qp](1)) + 
   0.6666666666666666667*_q44[_qp]*(-1.4142135623730950488*_u_x_grad[_qp](1) - 1.4142135623730950488*_u_y_grad[_qp](0) + _u_y_grad[_qp](2) + _u_z_grad[_qp](1)));
    }
    else if (jvar == _u_x_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](0)*_polar_z[_qp]*(13.333333333333334*_q11[_qp] + 26.666666666666668*_q12[_qp] - 26.666666666666668*_q44[_qp]) + 
   _grad_phi[_j][_qp](1)*_polar_y[_qp]*(9.428090415820634*_q11[_qp] - 9.428090415820634*_q12[_qp] - 18.856180831641268*_q44[_qp]) + 
   _grad_phi[_j][_qp](2)*_polar_x[_qp]*(13.333333333333334*_q11[_qp] - 13.333333333333334*_q12[_qp] + 13.333333333333334*_q44[_qp]) + 
   _grad_phi[_j][_qp](0)*_polar_x[_qp]*(-9.428090415820634*_q11[_qp] + 9.428090415820634*_q12[_qp] + 18.856180831641268*_q44[_qp]));
    }
    else if (jvar == _u_y_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](1)*_polar_z[_qp]*(13.333333333333334*_q11[_qp] + 26.666666666666668*_q12[_qp] - 26.666666666666668*_q44[_qp]) + 
   _grad_phi[_j][_qp](1)*_polar_x[_qp]*(9.428090415820634*_q11[_qp] - 9.428090415820634*_q12[_qp] - 18.856180831641268*_q44[_qp]) + 
   _grad_phi[_j][_qp](0)*_polar_y[_qp]*(9.428090415820634*_q11[_qp] - 9.428090415820634*_q12[_qp] - 18.856180831641268*_q44[_qp]) + 
   _grad_phi[_j][_qp](2)*_polar_y[_qp]*(13.333333333333334*_q11[_qp] - 13.333333333333334*_q12[_qp] + 13.333333333333334*_q44[_qp]));
    }
    else if (jvar == _u_z_var)
    {
      return -0.5*_test[_i][_qp] * (_grad_phi[_j][_qp](0)*_polar_x[_qp]*(13.333333333333334*_q11[_qp] - 13.333333333333334*_q12[_qp] + 13.333333333333334*_q44[_qp]) + 
   _grad_phi[_j][_qp](1)*_polar_y[_qp]*(13.333333333333334*_q11[_qp] - 13.333333333333334*_q12[_qp] + 13.333333333333334*_q44[_qp]) + 
   _grad_phi[_j][_qp](2)*_polar_z[_qp]*(13.333333333333334*_q11[_qp] + 26.666666666666668*_q12[_qp] + 53.333333333333336*_q44[_qp]));
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
