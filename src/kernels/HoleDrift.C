/**
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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/


#include "HoleDrift.h"

registerMooseObject("FerretApp", HoleDrift);

template<>
InputParameters validParams<HoleDrift>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a residual contribution due to nabla squared Phi = 0");
  params.addRequiredCoupledVar("p", "The variable p");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

HoleDrift::HoleDrift(const InputParameters & parameters)
  :Kernel(parameters),
    _p_var(coupled("p")),
    _p(coupledValue("p")),
   _len_scale(getParam<Real>("len_scale")),
   q(1.60217662 * 1e-19),
   kb(1.38064852 * 1e-23),
   T(300)
{
}

Real
HoleDrift::computeQpResidual()
{
  // Here u refers to the potential phi
  //  It is rescaled by q/kbT where
  //  q = 1.60217662 *1e-19
  //  kb = 1.38064852 *1e-23
  //  T = 300
  Real Relec = 0.0;
  Real k = q / (kb*T);
  Relec += -_p[_qp] *k* _grad_u[_qp] * _grad_test[_i][_qp] * _len_scale;
  ///  Moose::out << "\n R_elec-"; std::cout << " = " << Relec;
  return Relec;
}

Real
HoleDrift::computeQpJacobian()
{
  Real k = q / (kb*T);
  return -_grad_phi[_j][_qp] *_p[_qp]* k* _grad_test[_i][_qp] * _len_scale;
  // return 0.0;
}

Real
HoleDrift::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real k = q / (kb*T);
  if (jvar == _p_var)
    return -_grad_test[_i][_qp]* k* _phi[_j][_qp] * _grad_u[_qp] * _len_scale;
  else
    return 0.0;
}
