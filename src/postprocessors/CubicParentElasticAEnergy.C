/* AUTO-DERIVED from src/postprocessors/CubicParentElasticEnergy.C by substituting the polar order parameter for the
 * antiphase (AFD) tilt and the electrostrictive tensor Q for the rotostrictive
 * tensor R.  Same free energy, 1/2 (eps - eps0):C:(eps - eps0), with
 * eps0 = R.AA instead of Q.PP.  Generated 2026-09-14. */
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

#include "CubicParentElasticAEnergy.h"

registerMooseObject("FerretApp", CubicParentElasticAEnergy);

InputParameters CubicParentElasticAEnergy::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates a volume integral over the rotostrictive coupling energy density.");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase (AFD) tilt vector");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase (AFD) tilt vector");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the antiphase (AFD) tilt vector");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and J");
  return params;
}

CubicParentElasticAEnergy::CubicParentElasticAEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _antiphase_A_x(coupledValue("antiphase_A_x")),
   _antiphase_A_y(coupledValue("antiphase_A_y")),
   _antiphase_A_z(coupledValue("antiphase_A_z")),
   _R11(getMaterialProperty<Real>("R11")),
   _R12(getMaterialProperty<Real>("R12")),
   _R44(getMaterialProperty<Real>("R44")),
   _C11(getMaterialProperty<Real>("C11")),
   _C12(getMaterialProperty<Real>("C12")),
   _C44(getMaterialProperty<Real>("C44")),
   _energy_scale(getParam<Real>("energy_scale")),
   _strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_strain"))
{
}

Real
CubicParentElasticAEnergy::computeQpIntegral()
{
  return _energy_scale*((_C11[_qp]*Utility::pow<2>(Utility::pow<2>(_antiphase_A_x[_qp])*_R11[_qp] + (Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_R12[_qp] - _strain[_qp](0,0)) + 4*_C44[_qp]*Utility::pow<2>(-(_antiphase_A_x[_qp]*_antiphase_A_y[_qp]*_R44[_qp]) + _strain[_qp](0,1)) + 4*_C44[_qp]*Utility::pow<2>(-(_antiphase_A_x[_qp]*_antiphase_A_z[_qp]*_R44[_qp]) + _strain[_qp](0,2)) + 2*_C12[_qp]*(Utility::pow<2>(_antiphase_A_x[_qp])*_R11[_qp] + (Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_R12[_qp] - _strain[_qp](0,0))*(Utility::pow<2>(_antiphase_A_y[_qp])*_R11[_qp] + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_R12[_qp] - _strain[_qp](1,1)) + 
     _C11[_qp]*Utility::pow<2>(Utility::pow<2>(_antiphase_A_y[_qp])*_R11[_qp] + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_R12[_qp] - _strain[_qp](1,1)) + 4*_C44[_qp]*Utility::pow<2>(-(_antiphase_A_y[_qp]*_antiphase_A_z[_qp]*_R44[_qp]) + _strain[_qp](1,2)) + 2*_C12[_qp]*(Utility::pow<2>(_antiphase_A_x[_qp])*_R11[_qp] + (Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_R12[_qp] - _strain[_qp](0,0))*(Utility::pow<2>(_antiphase_A_z[_qp])*_R11[_qp] + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*_R12[_qp] - _strain[_qp](2,2)) + 2*_C12[_qp]*(Utility::pow<2>(_antiphase_A_y[_qp])*_R11[_qp] + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]))*_R12[_qp] - _strain[_qp](1,1))*(Utility::pow<2>(_antiphase_A_z[_qp])*_R11[_qp] + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*_R12[_qp] - _strain[_qp](2,2)) + 
     _C11[_qp]*Utility::pow<2>(Utility::pow<2>(_antiphase_A_z[_qp])*_R11[_qp] + (Utility::pow<2>(_antiphase_A_x[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp]))*_R12[_qp] - _strain[_qp](2,2)))/2.);
}
