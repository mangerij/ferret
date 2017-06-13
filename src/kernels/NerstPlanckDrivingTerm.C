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

#include "NerstPlanckDrivingTerm.h"

class NerstPlanckDrivingTerm;

template<>
InputParameters validParams<NerstPlanckDrivingTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("potential_int", "The electrostatic potential");
  params.addParam<Real>("mu_m", 1.0, "The mobility of the charge carriers");
  params.addParam<Real>("len_scale", 1.0, "The length scale of the unit");
  return params;
}

NerstPlanckDrivingTerm::NerstPlanckDrivingTerm(const InputParameters & parameters)
  :Kernel(parameters),
   _potential_int_var(coupled("potential_int")),
   _potential_int(coupledValue("potential_int")),
   _potential_int_grad(coupledGradient("potential_int")),
   _mu_m(getParam<Real>("mu_m")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
NerstPlanckDrivingTerm::computeQpResidual()
{
  // minus sign on the electric fields?
  return std::pow(_len_scale, 2.0) * _mu_m * (_u[_qp] * _potential_int_grad[_qp] * _grad_test[_i][_qp] + _grad_u[_qp] * _potential_int_grad[_qp] * _test[_i][_qp]);
}

Real
NerstPlanckDrivingTerm::computeQpJacobian()
{
  // minus sign on the electric fields?
  return std::pow(_len_scale, 2.0) * _mu_m * (_phi[_j][_qp] * _potential_int_grad[_qp] * _grad_test[_i][_qp] + _grad_phi[_j][_qp] * _potential_int_grad[_qp] * _test[_i][_qp]);
}

//need offdiag term for potential_int
