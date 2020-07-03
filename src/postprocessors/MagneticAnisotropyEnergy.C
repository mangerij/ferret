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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "MagneticAnisotropyEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MagneticAnisotropyEnergy);

template<>
InputParameters validParams<MagneticAnisotropyEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("potential_H_int", "The internal magnetic potential variable");
  params.addCoupledVar("potential_H_ext", 0.0, "The external magnetic potential variable");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the constrained magnetization");
  return params;
}

MagneticAnisotropyEnergy::MagneticAnisotropyEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _K1(getMaterialProperty<Real>("K1")),
   _K2(getMaterialProperty<Real>("K2")),
   _nx(getMaterialProperty<Real>("nx")),
   _ny(getMaterialProperty<Real>("ny")),
   _nz(getMaterialProperty<Real>("nz")),
   _Ms(getMaterialProperty<Real>("Ms"))
{
}

Real
MagneticAnisotropyEnergy::computeQpIntegral()
{
  return _K1[_qp]*Utility::pow<2>(_mag_x[_qp]*_nx[_qp] + _mag_y[_qp]*_ny[_qp] + _mag_z[_qp]*_nz[_qp])/(_Ms[_qp]*_Ms[_qp]);
}
