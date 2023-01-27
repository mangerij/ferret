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

#include "AFMSublatticeDMInteractionEnergyDensity.h"

registerMooseObject("FerretApp", AFMSublatticeDMInteractionEnergyDensity);

InputParameters AFMSublatticeDMInteractionEnergyDensity::validParams()
{

  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the DM interaction free energy density (coupling AFD and magnetic ordering).");
  params.addRequiredCoupledVar("mag1_x", "The x component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag1_y", "The y component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag1_z", "The z component of the constrained 1st sublattice magnetization vector");
  params.addRequiredCoupledVar("mag2_x", "The x component of the constrained 2nd sublattice magnetization vector");
  params.addRequiredCoupledVar("mag2_y", "The y component of the constrained 2nd sublattice magnetization vector");
  params.addRequiredCoupledVar("mag2_z", "The z component of the constrained 2nd sublattice magnetization vector");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase tilt vector");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase tilt vector");
  params.addRequiredCoupledVar("antiphase_A_z", "The z component of the antiphase tilt vector");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

AFMSublatticeDMInteractionEnergyDensity::AFMSublatticeDMInteractionEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
   _mag1_x(coupledValue("mag1_x")),
   _mag1_y(coupledValue("mag1_y")),
   _mag1_z(coupledValue("mag1_z")),
   _mag2_x(coupledValue("mag2_x")),
   _mag2_y(coupledValue("mag2_y")),
   _mag2_z(coupledValue("mag2_z")),
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
   _D0(getMaterialProperty<Real>("D0")),
   _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
AFMSublatticeDMInteractionEnergyDensity::computeValue()
{
  return  _energy_scale*(8.0*_D0[_qp]*(-(_antiphase_A_z[_qp]*_mag1_y[_qp]*_mag2_x[_qp]) + _antiphase_A_y[_qp]*_mag1_z[_qp]*_mag2_x[_qp] + _antiphase_A_z[_qp]*_mag1_x[_qp]*_mag2_y[_qp] - _antiphase_A_x[_qp]*_mag1_z[_qp]*_mag2_y[_qp] - _antiphase_A_y[_qp]*_mag1_x[_qp]*_mag2_z[_qp] + _antiphase_A_x[_qp]*_mag1_y[_qp]*_mag2_z[_qp]));
}
