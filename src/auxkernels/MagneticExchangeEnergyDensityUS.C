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

#include "MagneticExchangeEnergyDensityUS.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MagneticExchangeEnergyDensityUS);

InputParameters MagneticExchangeEnergyDensityUS::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("azimuthal_ph", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_th", "The polar component of the constrained magnetic vector");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and aJ");
  return params;
}

MagneticExchangeEnergyDensityUS::MagneticExchangeEnergyDensityUS(const InputParameters & parameters) :
  AuxKernel(parameters),
  _azimuthal_ph(coupledValue("azimuthal_ph")),
  _polar_th(coupledValue("polar_th")),
  _polar_th_grad(coupledGradient("polar_th")),
  _azimuthal_ph_grad(coupledGradient("azimuthal_ph")),
  _Ae(getMaterialProperty<Real>("Ae")),
  _Ms(getMaterialProperty<Real>("Ms")),
  _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
MagneticExchangeEnergyDensityUS::computeValue()
{
  return _energy_scale*((_Ae[_qp]*(Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)) + 2.0*(Utility::pow<2>(_polar_th_grad[_qp](0)) + Utility::pow<2>(_polar_th_grad[_qp](1)) + Utility::pow<2>(_polar_th_grad[_qp](2))) - (Utility::pow<2>(_azimuthal_ph_grad[_qp](0)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](1)) + Utility::pow<2>(_azimuthal_ph_grad[_qp](2)))*std::cos(2.0*_polar_th[_qp])))/2.);
}
