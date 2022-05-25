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

#include "AFMSingleIonCubicSixthAnisotropyEnergyDensity.h"

registerMooseObject("FerretApp", AFMSingleIonCubicSixthAnisotropyEnergyDensity);

InputParameters AFMSingleIonCubicSixthAnisotropyEnergyDensity::validParams()
{

  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates an integral over the DM interaction free energy density (coupling AFD and magnetic ordering).");
  params.addRequiredCoupledVar("mag_x", "The x component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag_y", "The y component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag_z", "The z component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase tilt vector");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase tilt vector");
  params.addRequiredCoupledVar("antiphase_A_z", "The z component of the antiphase tilt vector");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

AFMSingleIonCubicSixthAnisotropyEnergyDensity::AFMSingleIonCubicSixthAnisotropyEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
   _mag_x(coupledValue("mag_x")),
   _mag_y(coupledValue("mag_y")),
   _mag_z(coupledValue("mag_z")),
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
   _K1c(getMaterialProperty<Real>("K1c")),
   _Kt(getMaterialProperty<Real>("Kt")),
   _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
AFMSingleIonCubicSixthAnisotropyEnergyDensity::computeValue()
{
  return -_energy_scale*((_K1c[_qp] + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_Kt[_qp])*(Utility::pow<2>(_mag_y[_qp])*Utility::pow<2>(_mag_z[_qp])*Utility::pow<2>(_mag_x[_qp])));
}
