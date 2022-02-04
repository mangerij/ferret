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

#include "Transformed111RotostrictiveCouplingDispDerivative.h"
#include "libmesh/utility.h"

class Transformed111RotostrictiveCouplingDispDerivative;

registerMooseObject("FerretApp", Transformed111RotostrictiveCouplingDispDerivative);

InputParameters Transformed111RotostrictiveCouplingDispDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates a residual contribution due to the differentiation w.r.t spartial coordinates of the ferroelectric self-strain"
                             " in the condition for mechanical equilibrium. Note for BFO only.");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the tilt");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the tilt");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the tilt");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction in order parameter space this kernel acts in (e.g. for unrotated functionals 0 for q_x, 1 for q_y, 2 for q_z).");
  return params;
}

Transformed111RotostrictiveCouplingDispDerivative::Transformed111RotostrictiveCouplingDispDerivative(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _antiphase_A_x_var(coupled("antiphase_A_x")),
   _antiphase_A_y_var(coupled("antiphase_A_y")),
   _antiphase_A_z_var(coupled("antiphase_A_z")),
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
   _C11(getMaterialProperty<Real>("C11")),
   _C12(getMaterialProperty<Real>("C12")),
   _C44(getMaterialProperty<Real>("C44")),
   _R11(getMaterialProperty<Real>("R11")),
   _R12(getMaterialProperty<Real>("R12")),
   _R44(getMaterialProperty<Real>("R44"))
{
}

Real
Transformed111RotostrictiveCouplingDispDerivative::computeQpResidual()
{
  if (_component == 0)
  {
    return (0.16666666666666666667*(_C11[_qp]*((2.*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp]*(_antiphase_A_x[_qp] + 1.4142135623730950488*_antiphase_A_z[_qp]) + 
           _grad_test[_i][_qp](2)*(-1.4142135623730950488*Utility::pow<2>(_antiphase_A_x[_qp]) + 1.4142135623730950488*Utility::pow<2>(_antiphase_A_y[_qp]) + 4.*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]))*
         (_R11[_qp] - 1.*_R12[_qp]) + _grad_test[_i][_qp](0)*(2.8284271247461900976*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*(-1.*_R11[_qp] + _R12[_qp]) + 3.*Utility::pow<2>(_antiphase_A_x[_qp])*(_R11[_qp] + _R12[_qp]) + 
           2.*Utility::pow<2>(_antiphase_A_z[_qp])*(_R11[_qp] + 2.*_R12[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(_R11[_qp] + 5.*_R12[_qp]))) + 
     _C12[_qp]*((-1.4142135623730950488*_grad_test[_i][_qp](2)*Utility::pow<2>(_antiphase_A_y[_qp]) + _grad_test[_i][_qp](2)*_antiphase_A_x[_qp]*(1.4142135623730950488*_antiphase_A_x[_qp] - 4.*_antiphase_A_z[_qp]) - 
           2.*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp]*(_antiphase_A_x[_qp] + 1.4142135623730950488*_antiphase_A_z[_qp]))*(_R11[_qp] - 1.*_R12[_qp]) + 
        _grad_test[_i][_qp](0)*(2.8284271247461900976*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*(_R11[_qp] - 1.*_R12[_qp]) + 4.*Utility::pow<2>(_antiphase_A_z[_qp])*(_R11[_qp] + 2.*_R12[_qp]) + 
           3.*Utility::pow<2>(_antiphase_A_x[_qp])*(_R11[_qp] + 3.*_R12[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(5.*_R11[_qp] + 7.*_R12[_qp]))) + 
     4.*_C44[_qp]*(_grad_test[_i][_qp](1)*_antiphase_A_y[_qp]*(4.*_antiphase_A_x[_qp] - 2.8284271247461900976*_antiphase_A_z[_qp]) + 
        _grad_test[_i][_qp](2)*(1.4142135623730950488*Utility::pow<2>(_antiphase_A_x[_qp]) - 1.4142135623730950488*Utility::pow<2>(_antiphase_A_y[_qp]) + 2.*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]) + 
        _grad_test[_i][_qp](0)*(3.*Utility::pow<2>(_antiphase_A_x[_qp]) - 1.*Utility::pow<2>(_antiphase_A_y[_qp]) + 2.8284271247461900976*_antiphase_A_x[_qp]*_antiphase_A_z[_qp] - 2.*Utility::pow<2>(_antiphase_A_z[_qp])))*_R44[_qp]));
  }
  else if (_component == 1)
  {
    return (0.5*_grad_test[_i][_qp](1)*Utility::pow<2>(_antiphase_A_y[_qp])*(_C11[_qp]*(_R11[_qp] + _R12[_qp]) + _C12[_qp]*(_R11[_qp] + 3.*_R12[_qp]) + 4.*_C44[_qp]*_R44[_qp]) + 
   0.33333333333333333333*_antiphase_A_y[_qp]*(_C11[_qp]*(_grad_test[_i][_qp](0)*_antiphase_A_x[_qp] + 1.4142135623730950488*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] + 
         1.4142135623730950488*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp] + 2.*_grad_test[_i][_qp](2)*_antiphase_A_z[_qp])*(_R11[_qp] - 1.*_R12[_qp]) - 
      1.*_C12[_qp]*(_grad_test[_i][_qp](0)*_antiphase_A_x[_qp] + 1.4142135623730950488*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] + 1.4142135623730950488*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp] + 2.*_grad_test[_i][_qp](2)*_antiphase_A_z[_qp])*
       (_R11[_qp] - 1.*_R12[_qp]) + 4.*_C44[_qp]*(2.*_grad_test[_i][_qp](0)*_antiphase_A_x[_qp] - 1.4142135623730950488*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] - 
         1.4142135623730950488*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp] + _grad_test[_i][_qp](2)*_antiphase_A_z[_qp])*_R44[_qp]) + 
   0.16666666666666666667*_grad_test[_i][_qp](1)*(2.8284271247461900976*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*((_C11[_qp] - 1.*_C12[_qp])*(_R11[_qp] - 1.*_R12[_qp]) - 4.*_C44[_qp]*_R44[_qp]) + 
      2.*Utility::pow<2>(_antiphase_A_z[_qp])*((_C11[_qp] + 2.*_C12[_qp])*(_R11[_qp] + 2.*_R12[_qp]) - 4.*_C44[_qp]*_R44[_qp]) + 
      Utility::pow<2>(_antiphase_A_x[_qp])*(5.*_C12[_qp]*_R11[_qp] + 7.*_C12[_qp]*_R12[_qp] + _C11[_qp]*(_R11[_qp] + 5.*_R12[_qp]) - 4.*_C44[_qp]*_R44[_qp])));
  }
  else if (_component == 2)
  {
    return (0.16666666666666666667*(_C11[_qp]*(2.828427124746190098*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp]*(1.*_antiphase_A_x[_qp] + 1.414213562373095049*_antiphase_A_z[_qp])*(_R11[_qp] - 1.*_R12[_qp]) - 
        1.414213562373095049*_grad_test[_i][_qp](0)*(1.*Utility::pow<2>(_antiphase_A_x[_qp]) - 1.*Utility::pow<2>(_antiphase_A_y[_qp]) - 2.828427124746190098*_antiphase_A_x[_qp]*_antiphase_A_z[_qp])*
         (_R11[_qp] - 1.*_R12[_qp]) + 2.*_grad_test[_i][_qp](2)*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*(_R11[_qp] + 2.*_R12[_qp])) + 
     _C12[_qp]*(-2.828427124746190098*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp]*(1.*_antiphase_A_x[_qp] + 1.414213562373095049*_antiphase_A_z[_qp])*(_R11[_qp] - 1.*_R12[_qp]) + 
        _grad_test[_i][_qp](0)*(1.4142135623730950488*Utility::pow<2>(_antiphase_A_x[_qp]) - 1.4142135623730950488*Utility::pow<2>(_antiphase_A_y[_qp]) - 4.*_antiphase_A_x[_qp]*_antiphase_A_z[_qp])*(_R11[_qp] - 1.*_R12[_qp]) + 
        4.*_grad_test[_i][_qp](2)*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*(_R11[_qp] + 2.*_R12[_qp])) + 
     4.*_C44[_qp]*(_grad_test[_i][_qp](0)*(1.4142135623730950488*Utility::pow<2>(_antiphase_A_x[_qp]) - 1.4142135623730950488*Utility::pow<2>(_antiphase_A_y[_qp]) + 2.*_antiphase_A_x[_qp]*_antiphase_A_z[_qp]) - 
        2.*(_grad_test[_i][_qp](1)*_antiphase_A_y[_qp]*(1.4142135623730950488*_antiphase_A_x[_qp] - 1.*_antiphase_A_z[_qp]) + _grad_test[_i][_qp](2)*(Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]) - 2.*Utility::pow<2>(_antiphase_A_z[_qp]))))*
      _R44[_qp]));
  }
  else
    return 0.0;
}

Real
Transformed111RotostrictiveCouplingDispDerivative::computeQpJacobian()
{
  return 0.0;
}

Real
Transformed111RotostrictiveCouplingDispDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_component == 0)
  {
    if (jvar == _antiphase_A_x_var)
    {
      return _phi[_j][_qp] * (0.16666666666666666667*(_C11[_qp]*((-2.8284271247461900976*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] + 2.*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp] + 4.*_grad_test[_i][_qp](2)*_antiphase_A_z[_qp])*(_R11[_qp] - 1.*_R12[_qp]) + 
        2.8284271247461900976*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp]*(-1.*_R11[_qp] + _R12[_qp]) + 6.*_grad_test[_i][_qp](0)*_antiphase_A_x[_qp]*(_R11[_qp] + _R12[_qp])) + 
     _C12[_qp]*(2.8284271247461900976*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp]*(_R11[_qp] - 1.*_R12[_qp]) + 
        2.8284271247461901*(1.*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] - 0.707106781186547524*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp] - 1.41421356237309505*_grad_test[_i][_qp](2)*_antiphase_A_z[_qp])*
         (_R11[_qp] - 1.*_R12[_qp]) + 6.*_grad_test[_i][_qp](0)*_antiphase_A_x[_qp]*(_R11[_qp] + 3.*_R12[_qp])) + 
     24.*_C44[_qp]*(1.*_grad_test[_i][_qp](0)*_antiphase_A_x[_qp] + 0.47140452079103168*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] + 0.66666666666666667*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp] + 
        0.47140452079103168*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp] + 0.333333333333333333*_grad_test[_i][_qp](2)*_antiphase_A_z[_qp])*_R44[_qp]));
    }
    else if (jvar == _antiphase_A_y_var)
    {
      return _phi[_j][_qp] * (0.16666666666666666667*(_C12[_qp]*(-2.8284271247461900976*_grad_test[_i][_qp](2)*_antiphase_A_y[_qp] - 2.*_grad_test[_i][_qp](1)*(_antiphase_A_x[_qp] + 1.4142135623730950488*_antiphase_A_z[_qp]))*
      (_R11[_qp] - 1.*_R12[_qp]) + _C11[_qp]*(2.8284271247461900976*_grad_test[_i][_qp](2)*_antiphase_A_y[_qp] + 2.*_grad_test[_i][_qp](1)*(_antiphase_A_x[_qp] + 1.4142135623730950488*_antiphase_A_z[_qp]))*
      (_R11[_qp] - 1.*_R12[_qp]) + 10.*_C12[_qp]*_grad_test[_i][_qp](0)*_antiphase_A_y[_qp]*(1.*_R11[_qp] + 1.4*_R12[_qp]) + 2.*_C11[_qp]*_grad_test[_i][_qp](0)*_antiphase_A_y[_qp]*(_R11[_qp] + 5.*_R12[_qp]) + 
     16.*_C44[_qp]*(1.*_grad_test[_i][_qp](1)*_antiphase_A_x[_qp] - 0.5*_grad_test[_i][_qp](0)*_antiphase_A_y[_qp] - 0.707106781186547524*_grad_test[_i][_qp](2)*_antiphase_A_y[_qp] - 0.707106781186547524*_grad_test[_i][_qp](1)*_antiphase_A_z[_qp])*
      _R44[_qp]));
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return _phi[_j][_qp] * (0.16666666666666666667*(_C11[_qp]*((4.*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] + 2.8284271247461900976*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp])*(_R11[_qp] - 1.*_R12[_qp]) + 
        2.8284271247461900976*_grad_test[_i][_qp](0)*_antiphase_A_x[_qp]*(-1.*_R11[_qp] + _R12[_qp]) + 4.*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp]*(_R11[_qp] + 2.*_R12[_qp])) + 
     _C12[_qp]*(2.8284271247461900976*_grad_test[_i][_qp](0)*_antiphase_A_x[_qp]*(_R11[_qp] - 1.*_R12[_qp]) + 
        (-4.*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] - 2.8284271247461900976*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp])*(_R11[_qp] - 1.*_R12[_qp]) + 8.*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp]*(_R11[_qp] + 2.*_R12[_qp])) + 
     11.3137084989847604*_C44[_qp]*(1.*_grad_test[_i][_qp](0)*_antiphase_A_x[_qp] + 0.70710678118654752*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] - 1.*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp] - 
        1.41421356237309505*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp])*_R44[_qp]));
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
      return _phi[_j][_qp] * (0.33333333333333333333*_C11[_qp]*(_grad_test[_i][_qp](0) + 1.4142135623730950488*_grad_test[_i][_qp](2))*_antiphase_A_y[_qp]*(_R11[_qp] - 1.*_R12[_qp]) - 
   0.3333333333333333333*_C12[_qp]*(_grad_test[_i][_qp](0) + 1.4142135623730950488*_grad_test[_i][_qp](2))*_antiphase_A_y[_qp]*(_R11[_qp] - 1.*_R12[_qp]) + 
   2.666666666666666667*_C44[_qp]*(1.*_grad_test[_i][_qp](0) - 0.707106781186547524*_grad_test[_i][_qp](2))*_antiphase_A_y[_qp]*_R44[_qp] + 
   0.47140452079103168293*_grad_test[_i][_qp](1)*_antiphase_A_z[_qp]*((_C11[_qp] - 1.*_C12[_qp])*(_R11[_qp] - 1.*_R12[_qp]) - 4.*_C44[_qp]*_R44[_qp]) + 
   0.33333333333333333333*_grad_test[_i][_qp](1)*_antiphase_A_x[_qp]*(5.*_C12[_qp]*_R11[_qp] + 7.*_C12[_qp]*_R12[_qp] + _C11[_qp]*(_R11[_qp] + 5.*_R12[_qp]) - 4.*_C44[_qp]*_R44[_qp]));
    }
    else if (jvar == _antiphase_A_y_var)
    {
      return _phi[_j][_qp] * (0.33333333333333333333*_C11[_qp]*(_grad_test[_i][_qp](0)*_antiphase_A_x[_qp] + 1.4142135623730950488*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] + 1.4142135623730950488*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp] + 
      2.*_grad_test[_i][_qp](2)*_antiphase_A_z[_qp])*(_R11[_qp] - 1.*_R12[_qp]) - 0.3333333333333333333*_C12[_qp]*
    (_grad_test[_i][_qp](0)*_antiphase_A_x[_qp] + 1.4142135623730950488*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] + 1.4142135623730950488*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp] + 2.*_grad_test[_i][_qp](2)*_antiphase_A_z[_qp])*
    (_R11[_qp] - 1.*_R12[_qp]) + 1.3333333333333333333*_C44[_qp]*
    (2.*_grad_test[_i][_qp](0)*_antiphase_A_x[_qp] - 1.4142135623730950488*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] - 1.4142135623730950488*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp] + _grad_test[_i][_qp](2)*_antiphase_A_z[_qp])*_R44[_qp] + 
   1.*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp]*(_C11[_qp]*(_R11[_qp] + _R12[_qp]) + _C12[_qp]*(_R11[_qp] + 3.*_R12[_qp]) + 4.*_C44[_qp]*_R44[_qp]));
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return _phi[_j][_qp] * (0.471404520791031683*_C11[_qp]*(1.*_grad_test[_i][_qp](0) + 1.414213562373095049*_grad_test[_i][_qp](2))*_antiphase_A_y[_qp]*(_R11[_qp] - 1.*_R12[_qp]) - 
   0.471404520791031683*_C12[_qp]*(1.*_grad_test[_i][_qp](0) + 1.414213562373095049*_grad_test[_i][_qp](2))*_antiphase_A_y[_qp]*(_R11[_qp] - 1.*_R12[_qp]) + 
   1.3333333333333333333*_C44[_qp]*(-1.4142135623730950488*_grad_test[_i][_qp](0) + _grad_test[_i][_qp](2))*_antiphase_A_y[_qp]*_R44[_qp] + 
   0.47140452079103168293*_grad_test[_i][_qp](1)*_antiphase_A_x[_qp]*((_C11[_qp] - 1.*_C12[_qp])*(_R11[_qp] - 1.*_R12[_qp]) - 4.*_C44[_qp]*_R44[_qp]) + 
   0.6666666666666666667*_grad_test[_i][_qp](1)*_antiphase_A_z[_qp]*((_C11[_qp] + 2.*_C12[_qp])*(_R11[_qp] + 2.*_R12[_qp]) - 4.*_C44[_qp]*_R44[_qp]));
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
      return _phi[_j][_qp] * (0.16666666666666666667*(_C11[_qp]*(2.8284271247461900976*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp]*(_R11[_qp] - 1.*_R12[_qp]) - 
        2.828427124746190098*_grad_test[_i][_qp](0)*(1.*_antiphase_A_x[_qp] - 1.414213562373095049*_antiphase_A_z[_qp])*(_R11[_qp] - 1.*_R12[_qp]) + 4.*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp]*(_R11[_qp] + 2.*_R12[_qp])
        ) + _C12[_qp]*(-2.8284271247461900976*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp]*(_R11[_qp] - 1.*_R12[_qp]) + 
        _grad_test[_i][_qp](0)*(2.8284271247461900976*_antiphase_A_x[_qp] - 4.*_antiphase_A_z[_qp])*(_R11[_qp] - 1.*_R12[_qp]) + 8.*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp]*(_R11[_qp] + 2.*_R12[_qp])) + 
     11.3137084989847604*_C44[_qp]*(1.*_grad_test[_i][_qp](0)*_antiphase_A_x[_qp] - 1.41421356237309505*_grad_test[_i][_qp](2)*_antiphase_A_x[_qp] - 1.*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp] + 
        0.70710678118654752*_grad_test[_i][_qp](0)*_antiphase_A_z[_qp])*_R44[_qp]));
    }
    else if (jvar == _antiphase_A_y_var)
    {
      return _phi[_j][_qp] * (0.47140452079103168293*_C11[_qp]*_grad_test[_i][_qp](0)*_antiphase_A_y[_qp]*(_R11[_qp] - 1.*_R12[_qp]) - 0.47140452079103168293*_C12[_qp]*_grad_test[_i][_qp](0)*_antiphase_A_y[_qp]*(_R11[_qp] - 1.*_R12[_qp]) + 
   0.4714045207910316829*_C11[_qp]*_grad_test[_i][_qp](1)*(1.*_antiphase_A_x[_qp] + 1.414213562373095049*_antiphase_A_z[_qp])*(_R11[_qp] - 1.*_R12[_qp]) - 
   0.4714045207910316829*_C12[_qp]*_grad_test[_i][_qp](1)*(1.*_antiphase_A_x[_qp] + 1.414213562373095049*_antiphase_A_z[_qp])*(_R11[_qp] - 1.*_R12[_qp]) + 
   0.6666666666666666667*_C11[_qp]*_grad_test[_i][_qp](2)*_antiphase_A_y[_qp]*(_R11[_qp] + 2.*_R12[_qp]) + 1.3333333333333333333*_C12[_qp]*_grad_test[_i][_qp](2)*_antiphase_A_y[_qp]*(_R11[_qp] + 2.*_R12[_qp]) - 
   1.8856180831641267317*_C44[_qp]*_grad_test[_i][_qp](0)*_antiphase_A_y[_qp]*_R44[_qp] - 1.88561808316412673*_C44[_qp]*
    (1.*_grad_test[_i][_qp](1)*_antiphase_A_x[_qp] + 1.41421356237309505*_grad_test[_i][_qp](2)*_antiphase_A_y[_qp] - 0.70710678118654752*_grad_test[_i][_qp](1)*_antiphase_A_z[_qp])*_R44[_qp]);
    }
    else if (jvar == _antiphase_A_z_var)
    {
      return _phi[_j][_qp] * (0.6666666666666666667*_C11[_qp]*_grad_test[_i][_qp](0)*_antiphase_A_x[_qp]*(_R11[_qp] - 1.*_R12[_qp]) - 0.6666666666666666667*_C12[_qp]*_grad_test[_i][_qp](0)*_antiphase_A_x[_qp]*(_R11[_qp] - 1.*_R12[_qp]) + 
   0.6666666666666666667*_C11[_qp]*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp]*(_R11[_qp] - 1.*_R12[_qp]) - 0.6666666666666666667*_C12[_qp]*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp]*(_R11[_qp] - 1.*_R12[_qp]) + 
   0.6666666666666666667*_C11[_qp]*_grad_test[_i][_qp](2)*_antiphase_A_z[_qp]*(_R11[_qp] + 2.*_R12[_qp]) + 1.3333333333333333333*_C12[_qp]*_grad_test[_i][_qp](2)*_antiphase_A_z[_qp]*(_R11[_qp] + 2.*_R12[_qp]) + 
   1.3333333333333333333*_C44[_qp]*_grad_test[_i][_qp](0)*_antiphase_A_x[_qp]*_R44[_qp] + _C44[_qp]*
    (1.3333333333333333333*_grad_test[_i][_qp](1)*_antiphase_A_y[_qp] + 5.333333333333333333*_grad_test[_i][_qp](2)*_antiphase_A_z[_qp])*_R44[_qp]);
    }
    else
    {
      return 0.0;
    }
  }
  else
    return 0.0;
}
