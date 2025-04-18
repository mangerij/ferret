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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "MagnetostaticEnergyCart.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MagnetostaticEnergyCart);

InputParameters MagnetostaticEnergyCart::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("potential_H_int", "The internal magnetic potential variable");
  params.addCoupledVar("potential_H_ext", 0.0, "The external magnetic potential variable");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the constrained magnetization");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  params.addParam<Real>("mu0", 1.0, "permeability of the vacuum");
  return params;
}

MagnetostaticEnergyCart::MagnetostaticEnergyCart(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _potential_H_int_grad(coupledGradient("potential_H_int")),
   _potential_H_ext_grad(coupledGradient("potential_H_ext")),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _Ms(getMaterialProperty<Real>("Ms")),
   _mu0(getParam<Real>("mu0")),
   _energy_scale(getParam<Real>("energy_scale"))
{
  std::cout<<"__________________________________________________________________________"<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<" Magnetostatic Poisson equation:                                          "<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"       ∇·(∇ΦH)  = Ms ∇·m                                              "<<"\n";
  std::cout<<"__________________________________________________________________________"<<"\n";
}

Real
MagnetostaticEnergyCart::computeQpIntegral()
{
  // -1/2 * M*B = - 1/2 * M*(-gradPotential)
  return -_energy_scale*(0.5*_mu0*_Ms[_qp] * (-_potential_H_int_grad[_qp](0)*_mag_x[_qp]-_potential_H_int_grad[_qp](1)*_mag_y[_qp]-_potential_H_int_grad[_qp](2)*_mag_z[_qp]));
}
