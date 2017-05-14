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

#include "KappaTDiffusion.h"

template<>
InputParameters validParams<KappaTDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("temperature", "The local temperature");
  params.addParam<Real>("c0", 0.0, "the first fit coefficient");
  params.addParam<Real>("c1", 0.0, "the second fit coefficient");
  params.addParam<Real>("c2", 0.0, "the third fit coefficient");
  params.addParam<Real>("c3", 0.0, "the fourth fit coefficient");
  params.addParam<Real>("c4", 0.0, "the fifth fit coefficient");
//  params.addParam<Real>("c5", 0.0, "the sixth fit coefficient");
  return params;
}

KappaTDiffusion::KappaTDiffusion(const InputParameters & parameters)
  :Kernel(parameters),
   _temperature(coupledValue("temperature")),
   _temperature_grad(coupledGradient("temperature")),
   _c0(getParam<Real>("c0")),
   _c1(getParam<Real>("c1")),
   _c2(getParam<Real>("c2")),
   _c3(getParam<Real>("c3")),
   _c4(getParam<Real>("c4"))
{
}

Real
KappaTDiffusion::computeQpResidual()
{
  return (_c0 * _temperature[_qp] *_temperature[_qp] * std::exp(-_c1 *_temperature[_qp])+_c2*_temperature[_qp]*std::exp(- _c3 *_temperature[_qp])  + _c4 * _temperature[_qp] ) * _temperature_grad[_qp] * _grad_test[_i][_qp];
}

Real
KappaTDiffusion::computeQpJacobian()
{
  return (_c4 + _c2 *std::exp(-_c3 *_temperature[_qp]) + 2.0 * _c0 * std::exp(-_c1 *_temperature[_qp]) *_temperature[_qp] - _c2 *_c3 *std::exp(-_c3 *_temperature[_qp]) *_temperature[_qp] -  _c0* _c1 *std::exp(-_c1 *_temperature[_qp])* _temperature[_qp] *_temperature[_qp])*_phi[_j][_qp] *_temperature_grad[_qp]* _grad_test[_i][_qp] + _grad_phi[_j][_qp]* _grad_test[_i][_qp] *(_c0 * _temperature[_qp] *_temperature[_qp] * std::exp(-_c1 *_temperature[_qp])+_c2*_temperature[_qp]*std::exp(- _c3 *_temperature[_qp])  + _c4 * _temperature[_qp]);
}
