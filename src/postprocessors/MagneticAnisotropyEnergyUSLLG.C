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

   For help with FERRET please contact J. Mangeri <john.m.mangeri@gmail.com>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "MagneticAnisotropyEnergyUSLLG.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MagneticAnisotropyEnergyUSLLG);

InputParameters MagneticAnisotropyEnergyUSLLG::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over the magnetic exchange energy density.");
  params.addRequiredCoupledVar("azimuthal_ph", "The azimuthal component of the constrained magnetic vector");
  params.addRequiredCoupledVar("polar_th", "The polar component of the constrained magnetic vector");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and aJ");
  return params;
}

MagneticAnisotropyEnergyUSLLG::MagneticAnisotropyEnergyUSLLG(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  //  _component(getParam<unsigned int>("component")),
  _polar_th(coupledValue("polar_th")),
  _azimuthal_ph(coupledValue("azimuthal_ph")),
  _K1(getMaterialProperty<Real>("K1")),
   //_K2(getMaterialProperty<Real>("K2")),
  _nx(getMaterialProperty<Real>("nx")),
  _ny(getMaterialProperty<Real>("ny")),
  _nz(getMaterialProperty<Real>("nz")),
  _Ms(getMaterialProperty<Real>("Ms")),
  _mu0(getMaterialProperty<Real>("mu0"))
{
}

Real
MagneticAnisotropyEnergyUSLLG::computeQpIntegral()
{
  return -_K1[_qp]*Utility::pow<2>(std::cos( _azimuthal_ph[_qp])*std::sin(_polar_th[_qp])*_nx[_qp]
				            +std::sin( _azimuthal_ph[_qp])*std::sin(_polar_th[_qp])*_ny[_qp] + std::cos(_polar_th[_qp])*_nz[_qp]);
}
