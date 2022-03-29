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

#include "ElectrostaticEnergy.h"

registerMooseObject("FerretApp", ElectrostaticEnergy);

InputParameters ElectrostaticEnergy::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over the P$*$E term.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("potential_E_int", "The internal electric potential");
  params.addCoupledVar("potential_E_ext", 0.0, "The external electric potential");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

ElectrostaticEnergy::ElectrostaticEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _potential_E_int_grad(coupledGradient("potential_E_int")),
  _potential_E_ext_grad(coupledGradient("potential_E_ext")),
  _len_scale(getParam<Real>("len_scale")),
  _energy_scale(getParam<Real>("energy_scale"))
{
  std::cout<<"__________________________________________________________________________"<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<" Electrostatic Poisson equation:                                          "<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"       ∇·(ε∇Φ)  = -∇·P                                                   "<<"\n";
  std::cout<<"__________________________________________________________________________"<<"\n";
}

Real
ElectrostaticEnergy::computeQpIntegral()
{
  RealVectorValue P;
  P(0) = _polar_x[_qp]; P(1) = _polar_y[_qp]; P(2) = _polar_z[_qp];
  return _energy_scale*((0.5 * P * _potential_E_int_grad[_qp]) * std::pow(_len_scale, 2.0) + (P * _potential_E_ext_grad[_qp]) * std::pow(_len_scale, 2.0));
}
