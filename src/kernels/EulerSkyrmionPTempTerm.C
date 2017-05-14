/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#include "EulerSkyrmionPTempTerm.h"

class EulerSkyrmionPTempTerm;

template<>
InputParameters validParams<EulerSkyrmionPTempTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("theta", "The theta variable");
  params.addRequiredCoupledVar("P", "The polar magnitude variable");
  params.addParam<Real>("t", 1.0, "the alpha1 temperature constant");
  params.addParam<Real>("kappa", 1.0, "the anisotropy constant");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("xi0", 1.0, "the domain wall coefficient");
  return params;
}

EulerSkyrmionPTempTerm::EulerSkyrmionPTempTerm(const InputParameters & parameters)
  :Kernel(parameters),
   _theta_var(coupled("theta")),
   _P_var(coupled("P")),
   _theta(coupledValue("theta")),
   _theta_grad(coupledGradient("theta")),
   _P(coupledValue("P")),
   _t(getParam<Real>("t")),
   _kappa(getParam<Real>("kappa")),
   _len_scale(getParam<Real>("len_scale")),
   _xi0(getParam<Real>("xi0"))
{
}

Real
EulerSkyrmionPTempTerm::computeQpResidual()
{
  return _test[_i][_qp] * (
  _t + _xi0 * _xi0 * _theta_grad[_qp](0) * _theta_grad[_qp](0)
  + (_kappa + _xi0 * _xi0 / (_q_point[_qp](0) * _q_point[_qp](0))) * std::sin(_theta[_qp]) * std::sin(_theta[_qp]) ) * _P[_qp];
}

Real
EulerSkyrmionPTempTerm::computeQpJacobian()
{
  return _test[_i][_qp] * (
  _t + _xi0 * _xi0 * _theta_grad[_qp](0) * _theta_grad[_qp](0)
  + (_kappa + _xi0 * _xi0 / (_q_point[_qp](0) * _q_point[_qp](0))) * std::sin(_theta[_qp]) * std::sin(_theta[_qp]) ) * _phi[_j][_qp];
}

Real
EulerSkyrmionPTempTerm::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
