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

#include "AFMSublatticeAnisotropyEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", AFMSublatticeAnisotropyEnergy);

InputParameters AFMSublatticeAnisotropyEnergy::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the constrained magnetization");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

AFMSublatticeAnisotropyEnergy::AFMSublatticeAnisotropyEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _K1(getMaterialProperty<Real>("K1")),
   _Ms(getMaterialProperty<Real>("Ms")),
  _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
AFMSublatticeAnisotropyEnergy::computeQpIntegral()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  RealVectorValue f = w/std::sqrt(w*w);
  return _energy_scale*(-_K1[_qp]*Utility::pow<2>(_mag_x[_qp]*f(0) + _mag_y[_qp]*f(1) + _mag_z[_qp]*f(2)));//*(_Ms[_qp]*_Ms[_qp]);
}
