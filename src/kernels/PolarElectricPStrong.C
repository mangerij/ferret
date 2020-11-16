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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "PolarElectricPStrong.h"

class PolarElectricPStrong;

registerMooseObject("FerretApp", PolarElectricPStrong);

template<>
InputParameters validParams<PolarElectricPStrong>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to -P$*$E term in the total energy.");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  params.addRequiredCoupledVar("potential_E_int", "The internal electrostatic potential variable");
  params.addCoupledVar("potential_E_ext", 0.0, "The external electrostatic potential variable");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

PolarElectricPStrong::PolarElectricPStrong(const InputParameters & parameters)
  :Kernel(parameters),
   _component(getParam<unsigned int>("component")),
   _potential_E_int_var(coupled("potential_E_int")),
   _potential_E_ext_var(coupled("potential_E_ext")),
   _potential_E_int_grad(coupledGradient("potential_E_int")),
   _potential_E_ext_grad(coupledGradient("potential_E_ext")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
PolarElectricPStrong::computeQpResidual()
{
    Real RpolarP = 0.0;

    RpolarP += (_potential_E_int_grad[_qp](_component) + _potential_E_ext_grad[_qp](_component)) * _test[_i][_qp] * std::pow(_len_scale, 2.0);

    ///  Moose::out << "\n R_polarP-"; std::cout << _component << " = " << RpolarP;

    return RpolarP;
}

Real
PolarElectricPStrong::computeQpJacobian()
{
  return 0.0;
}

Real
PolarElectricPStrong::computeQpOffDiagJacobian(unsigned int jvar)
{
    if( jvar == _potential_E_int_var )
      return  _grad_phi[_j][_qp](_component) * _test[_i][_qp] * std::pow(_len_scale, 2.0);
    else if( jvar == _potential_E_ext_var)
      return  _grad_phi[_j][_qp](_component) * _test[_i][_qp] * std::pow(_len_scale, 2.0);
    else
    {
      return 0.0;
    }
}
