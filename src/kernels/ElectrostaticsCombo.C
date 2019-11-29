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

#include "ElectrostaticsCombo.h"

registerMooseObject("FerretApp", ElectrostaticsCombo);

template<>
InputParameters validParams<ElectrostaticsCombo>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  params.addRequiredParam<Real>("permittivity", "permittivity");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  params.addRequiredCoupledVar("n", "The variable n");
  params.addRequiredCoupledVar("p", "The variable p");
  return params;
}

ElectrostaticsCombo::ElectrostaticsCombo(const InputParameters & parameters)
  :Kernel(parameters),
   _permittivity(getParam<Real>("permittivity")),
   _len_scale(getParam<Real>("len_scale")),
   _n(coupledValue("n")),
   _p(coupledValue("p")),
   _n_var(coupled("n")),
   _p_var(coupled("p")),
   q(1.60217662 * 1e-19),
   kb(1.38064852 * 1e-23),
   T(300)
{
}

Real
ElectrostaticsCombo::computeQpResidual()
{
  Real k = q / (kb*T);
  Real Relec = 0.0;
  Relec += _grad_test[_i][_qp] *k* (_permittivity * _grad_u[_qp] *  _len_scale);
  Relec += _test[_i][_qp] * (_p[_qp] - _n[_qp]);
  ///  Moose::out << "\n R_elec-"; std::cout << " = " << Relec;
  return Relec;
}

Real
ElectrostaticsCombo::computeQpJacobian()
{
  Real k = q / (kb*T);
  return _permittivity *k* _grad_phi[_j][_qp] * _grad_test[_i][_qp] * _len_scale;
}

Real
ElectrostaticsCombo::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _n_var)
    return -_test[_i][_qp]*_phi[_j][_qp];
  else if (jvar == _p_var)
    return _test[_i][_qp]*_phi[_j][_qp];
  else
    return 0.0;
}
