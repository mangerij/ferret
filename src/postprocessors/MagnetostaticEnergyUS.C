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

#include "MagnetostaticEnergyUS.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MagnetostaticEnergyUS);

InputParameters MagnetostaticEnergyUS::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("potential_H_int", "The internal magnetic potential variable");
  params.addCoupledVar("potential_H_ext", 0.0, "The external magnetic potential variable");
  params.addRequiredCoupledVar("azimuth_phi", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_theta", "The polar component of the constrained magnetic vector");
  return params;
}

MagnetostaticEnergyUS::MagnetostaticEnergyUS(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _potential_H_int_grad(coupledGradient("potential_H_int")),
   _potential_H_ext_grad(coupledGradient("potential_H_ext")),
   _azimuth_phi(coupledValue("azimuth_phi")),
   _polar_theta(coupledValue("polar_theta")),
   _Ms(getMaterialProperty<Real>("Ms"))
{
}

Real
MagnetostaticEnergyUS::computeQpIntegral()
{
  // -1/2 * M*B = - 1/2 * M*(-gradPotential)
  return -0.5*_Ms[_qp] * (-_potential_H_int_grad[_qp](0)*std::cos(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]) - _potential_H_int_grad[_qp](1)*std::sin(_azimuth_phi[_qp])*std::sin(_polar_theta[_qp]) - _potential_H_int_grad[_qp](2)*std::cos(_polar_theta[_qp]));
}
