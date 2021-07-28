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

#include "CorrectionFactorDistortBFOEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", CorrectionFactorDistortBFOEnergy);

template<>
InputParameters validParams<CorrectionFactorDistortBFOEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral whose integrand is the correction to the local free energy.");
  params.addRequiredCoupledVar("antiferrodis_A_x", "The x component of the antiferrodistortive tilt");
  params.addRequiredCoupledVar("antiferrodis_A_y", "The y component of the antiferrodistortive tilt");
  params.addCoupledVar("antiferrodis_A_z", 0.0, "The z component of the antiferrodistortive tilt");
  return params;
}

CorrectionFactorDistortBFOEnergy::CorrectionFactorDistortBFOEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _antiferrodis_A_x(coupledValue("antiferrodis_A_x")),
  _antiferrodis_A_y(coupledValue("antiferrodis_A_y")),
  _antiferrodis_A_z(coupledValue("antiferrodis_A_z")),
  _f1(getMaterialProperty<Real>("f1"))
{
}

Real
CorrectionFactorDistortBFOEnergy::computeQpIntegral()
{
  return (_f1[_qp]*(-Utility::pow<4>(_antiferrodis_A_x[_qp]) + _antiferrodis_A_x[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_y[_qp] - Utility::pow<4>(_antiferrodis_A_y[_qp]) + _antiferrodis_A_x[_qp]*_antiferrodis_A_x[_qp]*_antiferrodis_A_z[_qp]*_antiferrodis_A_z[_qp] + _antiferrodis_A_y[_qp]*_antiferrodis_A_y[_qp]*_antiferrodis_A_z[_qp]*_antiferrodis_A_z[_qp] - Utility::pow<4>(_antiferrodis_A_z[_qp])));
}
