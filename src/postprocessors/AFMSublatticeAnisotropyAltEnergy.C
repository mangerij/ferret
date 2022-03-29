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

#include "AFMSublatticeAnisotropyAltEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", AFMSublatticeAnisotropyAltEnergy);

InputParameters AFMSublatticeAnisotropyAltEnergy::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("mag1_x", "The x component of the constrained magnetization");
  params.addRequiredCoupledVar("mag1_y", "The y component of the constrained magnetization");
  params.addCoupledVar("mag1_z", 0.0, "The z component of the constrained magnetization");
  params.addCoupledVar("mag2_x", 0.0, "The x component of the constrained magnetization");
  params.addCoupledVar("mag2_y", 0.0, "The y component of the constrained magnetization");
  params.addCoupledVar("mag2_z", 0.0, "The z component of the constrained magnetization");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

AFMSublatticeAnisotropyAltEnergy::AFMSublatticeAnisotropyAltEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _mag1_x(coupledValue("mag1_x")),
   _mag1_y(coupledValue("mag1_y")),
   _mag1_z(coupledValue("mag1_z")),
   _mag2_x(coupledValue("mag2_x")),
   _mag2_y(coupledValue("mag2_y")),
   _mag2_z(coupledValue("mag2_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _K1(getMaterialProperty<Real>("K1")),
   _Ms(getMaterialProperty<Real>("Ms")),
  _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
AFMSublatticeAnisotropyAltEnergy::computeQpIntegral()
{
  RealVectorValue w(_polar_x[_qp], _polar_y[_qp], _polar_z[_qp]);
  RealVectorValue f = w/std::sqrt(w*w);
  return _energy_scale*(-(_K1[_qp]*Utility::pow<2>(_mag1_x[_qp]*f(0) + _mag1_y[_qp]*f(1) + _mag1_z[_qp]*f(2)))-(_K1[_qp]*Utility::pow<2>(_mag1_x[_qp]*f(0) + _mag1_y[_qp]*f(1) + _mag1_z[_qp]*f(2))));
}
