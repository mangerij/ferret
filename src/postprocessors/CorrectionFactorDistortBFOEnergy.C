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

#include "CorrectionFactorDistortBFOEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", CorrectionFactorDistortBFOEnergy);

InputParameters CorrectionFactorDistortBFOEnergy::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral whose integrand is the correction to the local free energy.");
  params.addRequiredCoupledVar("antiphase_A_x", "The x component of the antiphase tilt");
  params.addRequiredCoupledVar("antiphase_A_y", "The y component of the antiphase tilt");
  params.addCoupledVar("antiphase_A_z", 0.0, "The z component of the antiphase tilt");
  return params;
}

CorrectionFactorDistortBFOEnergy::CorrectionFactorDistortBFOEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _antiphase_A_x(coupledValue("antiphase_A_x")),
  _antiphase_A_y(coupledValue("antiphase_A_y")),
  _antiphase_A_z(coupledValue("antiphase_A_z")),
  _f1(getMaterialProperty<Real>("f1"))
{
}

Real
CorrectionFactorDistortBFOEnergy::computeQpIntegral()
{
  return ((-Utility::pow<6>(_antiphase_A_x[_qp]) - Utility::pow<6>(_antiphase_A_y[_qp]) + Utility::pow<4>(_antiphase_A_y[_qp])*(-1 + Utility::pow<2>(_antiphase_A_z[_qp])) - Utility::pow<4>(_antiphase_A_z[_qp])*(1 + Utility::pow<2>(_antiphase_A_z[_qp])) + Utility::pow<4>(_antiphase_A_x[_qp])*(-1 + Utility::pow<2>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp])) + 
     Utility::pow<2>(_antiphase_A_y[_qp])*(Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp])) + Utility::pow<2>(_antiphase_A_x[_qp])*(Utility::pow<4>(_antiphase_A_y[_qp]) + Utility::pow<2>(_antiphase_A_z[_qp]) + Utility::pow<4>(_antiphase_A_z[_qp]) + Utility::pow<2>(_antiphase_A_y[_qp])*(1 - 3*Utility::pow<2>(_antiphase_A_z[_qp]))))*_f1[_qp]);
}
