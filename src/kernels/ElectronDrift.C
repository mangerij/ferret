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

#include "ElectronDrift.h"

registerMooseObject("FerretApp", ElectronDrift);

template<>
InputParameters validParams<ElectronDrift>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  params.addRequiredCoupledVar("n", "The variable n or p");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

ElectronDrift::ElectronDrift(const InputParameters & parameters)
  :Kernel(parameters),
    _n_var(coupled("n")),
    _n(coupledValue("n")),
   _len_scale(getParam<Real>("len_scale")),
   q(1.60217662 * 1e-19),
   kb(1.38064852 * 1e-23),
   T(300)
{
}

Real
ElectronDrift::computeQpResidual()
{
  // Here u refers to the potential phi
  //  It is rescaled by q/kbT where
  //  q = 1.60217662 *1e-19
  //  kb = 1.38064852 *1e-23
  //  T = 300
  Real Relec = 0.0;
  Real k = q / (kb*T);
  Relec += _n[_qp] *k* _grad_u[_qp] * _grad_test[_i][_qp] * _len_scale;
  ///  Moose::out << "\n R_elec-"; std::cout << " = " << Relec;
  return Relec;
}

Real
ElectronDrift::computeQpJacobian()
{
  Real k = q / (kb*T);
  return _grad_phi[_j][_qp] *_n[_qp]* k* _grad_test[_i][_qp] * _len_scale;
}

Real
ElectronDrift::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real k = q / (kb*T);
  if (jvar == _n_var)
    return _grad_test[_i][_qp]* k* _phi[_j][_qp] * _grad_u[_qp] * _len_scale;
  else
    return 0.0;
}
