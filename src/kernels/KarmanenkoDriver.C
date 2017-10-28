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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "KarmanenkoDriver.h"

template<>
InputParameters validParams<KarmanenkoDriver>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
  params.addCoupledVar("potential_ext", 0.0, "The external electric potential variable");
  params.addRequiredCoupledVar("temperature", "The temperature at the grid point");
  params.addParam<Real>("rho1", 1.0, "the density of the electrocaloric material");
  params.addParam<Real>("C1", 1.0, "the first coefficient of the residual contribution");
  params.addParam<Real>("C2", 1.0, "the second coefficient of the residual contribution");
  params.addParam<Real>("C3", 1.0, "the third coefficient of the residual contribution");
  params.addParam<Real>("C4", 1.0, "the fourth coefficient of the residual contribution");
  params.addParam<Real>("dEstep", 0.0, "change in the electric field as a function of time");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

KarmanenkoDriver::KarmanenkoDriver(const InputParameters & parameters)
  :Kernel(parameters),
   _potential_int_var(coupled("potential_int")),
   _potential_ext_var(coupled("potential_ext")),
   _potential_int_grad(coupledGradient("potential_int")),
   _potential_ext_grad(coupledGradient("potential_ext")),
   _temperature_var(coupled("temperature")),
   _temperature(coupledValue("temperature")),
   _rho1(getParam<Real>("rho1")),
   _C1(getParam<Real>("C1")),
   _C2(getParam<Real>("C2")),
   _C3(getParam<Real>("C3")),
   _C4(getParam<Real>("C4")),
   _dEstep(getParam<Real>("dEstep")),
   _len_scale(getParam<Real>("len_scale"))
{
  std::cout<<"Implementing Karmanenko field-induced entropic change with dEstep = "<< _dEstep <<"\n";
}

Real
KarmanenkoDriver::computeQpResidual()

////0.000116379 Ez + 0.00117004 Ez^2 - 3.27924*10^-6 Ez T + 4.68114*10^-6 Ez^2 T
{
  return _test[_i][_qp] * std::pow(_len_scale, 2.0) * _rho1 * _temperature[_qp] * (- _C1 * _potential_int_grad[_qp](2) + _C2 * _potential_int_grad[_qp](2) * _potential_int_grad[_qp](2) - _C3 * _potential_int_grad[_qp](2) * _temperature[_qp] + _C4 * _potential_int_grad[_qp](2) * _potential_int_grad[_qp](2) * _temperature[_qp]) * _dEstep;
}

Real
KarmanenkoDriver::computeQpJacobian()
{
  return _test[_i][_qp] * std::pow(_len_scale, 2.0) * _rho1 * _phi[_j][_qp] * (- _C1 * _potential_int_grad[_qp](2) + _C2 * _potential_int_grad[_qp](2) * _potential_int_grad[_qp](2) - _C3 * _potential_int_grad[_qp](2) * _temperature[_qp] + _C4 * _potential_int_grad[_qp](2) * _potential_int_grad[_qp](2) * _temperature[_qp]) * _dEstep;
}

Real
KarmanenkoDriver::computeQpOffDiagJacobian(unsigned int jvar)
{
    if( jvar == _potential_int_var )
      return  _test[_i][_qp] * std::pow(_len_scale, 2.0) * _rho1 * _temperature[_qp] * (- _C1 * _grad_phi[_j][_qp](2) + _C2 * _grad_phi[_j][_qp](2) * _potential_int_grad[_qp](2) - _C3 * _grad_phi[_j][_qp](2) * _temperature[_qp] + _C4 * _grad_phi[_j][_qp](2) * _potential_int_grad[_qp](2) * _temperature[_qp]) * _dEstep;
    else if( jvar == _potential_ext_var)
      return  0.0;
    else
    {
      return 0.0;
    }
}
