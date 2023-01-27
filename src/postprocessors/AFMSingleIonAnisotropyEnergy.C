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

#include "AFMSingleIonAnisotropyEnergy.h"

registerMooseObject("FerretApp", AFMSingleIonAnisotropyEnergy);

InputParameters AFMSingleIonAnisotropyEnergy::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over the DM interaction free energy density (coupling AFD and magnetic ordering).");
  params.addRequiredCoupledVar("mag1_x", "The x component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag1_y", "The y component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag1_z", "The z component of the constrained 1st sublattice magnetization vector");
  params.addCoupledVar("mag2_x", 0.0, "The x component of the constrained 2nd sublattice magnetization vector");
  params.addCoupledVar("mag2_y", 0.0, "The y component of the constrained 2nd sublattice magnetization vector");
  params.addCoupledVar("mag2_z", 0.0, "The z component of the constrained 2nd sublattice magnetization vector");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

AFMSingleIonAnisotropyEnergy::AFMSingleIonAnisotropyEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _mag1_x(coupledValue("mag1_x")),
   _mag1_y(coupledValue("mag1_y")),
   _mag1_z(coupledValue("mag1_z")),
   _mag2_x(coupledValue("mag2_x")),
   _mag2_y(coupledValue("mag2_y")),
   _mag2_z(coupledValue("mag2_z")),
   _K1c(getMaterialProperty<Real>("K1c")),
   _K2c(getMaterialProperty<Real>("K2c")),
   _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
AFMSingleIonAnisotropyEnergy::computeQpIntegral()
{
  return _energy_scale*(_K2c[_qp]*Utility::pow<2>(_mag1_x[_qp] + _mag2_x[_qp])*Utility::pow<2>(_mag1_y[_qp] + _mag2_y[_qp])*Utility::pow<2>(_mag1_z[_qp] + _mag2_z[_qp]) + _K1c[_qp]*(Utility::pow<2>(_mag1_z[_qp])*Utility::pow<2>(_mag2_x[_qp]) + Utility::pow<2>(_mag1_z[_qp])*Utility::pow<2>(_mag2_y[_qp]) + Utility::pow<2>(_mag2_x[_qp])*Utility::pow<2>(_mag2_y[_qp]) + 2.0*_mag1_z[_qp]*(Utility::pow<2>(_mag2_x[_qp]) + Utility::pow<2>(_mag2_y[_qp]))*_mag2_z[_qp] + 
      (Utility::pow<2>(_mag2_x[_qp]) + Utility::pow<2>(_mag2_y[_qp]))*Utility::pow<2>(_mag2_z[_qp]) + Utility::pow<2>(_mag1_y[_qp])*(Utility::pow<2>(_mag2_x[_qp]) + Utility::pow<2>(_mag1_z[_qp] + _mag2_z[_qp])) + 2.0*_mag1_y[_qp]*_mag2_y[_qp]*(Utility::pow<2>(_mag2_x[_qp]) + Utility::pow<2>(_mag1_z[_qp] + _mag2_z[_qp])) + Utility::pow<2>(_mag1_x[_qp])*(Utility::pow<2>(_mag1_y[_qp] + _mag2_y[_qp]) + Utility::pow<2>(_mag1_z[_qp] + _mag2_z[_qp])) + 2.0*_mag1_x[_qp]*_mag2_x[_qp]*(Utility::pow<2>(_mag1_y[_qp] + _mag2_y[_qp]) + Utility::pow<2>(_mag1_z[_qp] + _mag2_z[_qp]))));
}
