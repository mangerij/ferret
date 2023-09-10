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

#include "MagneticExcessLLBEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MagneticExcessLLBEnergy);

InputParameters MagneticExcessLLBEnergy::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the constrained magnetization");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

MagneticExcessLLBEnergy::MagneticExcessLLBEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _alpha(getMaterialProperty<Real>("alpha")),
  _alpha_long(getMaterialProperty<Real>("alpha_long")),
  _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
MagneticExcessLLBEnergy::computeQpIntegral()
{
  return _energy_scale*((0.25)*(_alpha_long[_qp]*Utility::pow<2>(_mag_x[_qp]*_mag_x[_qp]+_mag_y[_qp]*_mag_y[_qp]+_mag_z[_qp]*_mag_z[_qp]-1.0))/(1.0+_alpha[_qp]*_alpha[_qp]));
}
