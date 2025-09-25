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

#include "MasterMagneticZeemanEnergyCart.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MasterMagneticZeemanEnergyCart);

InputParameters MasterMagneticZeemanEnergyCart::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
 params.addClassDescription("Calculates a volume integral over the Zeeman interaction energy.");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the constrained magnetization");
  params.addRequiredCoupledVar("Hext_x", "The fixed x-component of the magnetic field variable");
  params.addRequiredCoupledVar("Hext_y", "The fixed y-component of the magnetic field variable");
  params.addRequiredCoupledVar("Hext_z", "The fixed z-component of the magnetic field variable");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  params.addParam<Real>("Hscale", 1.0, "scale for magnetic field");
  params.addParam<Real>("mu0", 1.0, "permeability of the vacuum");
  return params;
}

MasterMagneticZeemanEnergyCart::MasterMagneticZeemanEnergyCart(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _Hext_x(coupledValue("Hext_x")),
   _Hext_y(coupledValue("Hext_y")),
   _Hext_z(coupledValue("Hext_z")),
   _energy_scale(getParam<Real>("energy_scale")),
   _Ms(getMaterialProperty<Real>("Ms")),
   _Hscale(getParam<Real>("Hscale")),
   _mu0(getParam<Real>("mu0"))
{
}

Real
MasterMagneticZeemanEnergyCart::computeQpIntegral()
{
  return -_energy_scale*(_mu0*_Ms[_qp] * (_Hext_x[_qp]*_mag_x[_qp]+_Hext_y[_qp]*_mag_y[_qp]+_Hext_z[_qp]*_mag_z[_qp]));
}
