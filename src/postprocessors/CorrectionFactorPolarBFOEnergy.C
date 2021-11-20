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

   You should have received a co_polar_y[_qp] of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "CorrectionFactorPolarBFOEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", CorrectionFactorPolarBFOEnergy);

InputParameters CorrectionFactorPolarBFOEnergy::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral whose integrand is the correction to the local free energy.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

CorrectionFactorPolarBFOEnergy::CorrectionFactorPolarBFOEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _f0(getMaterialProperty<Real>("f0"))
{
}

Real
CorrectionFactorPolarBFOEnergy::computeQpIntegral()
{
  return (_f0[_qp]*(-Utility::pow<6>(_polar_x[_qp]) - Utility::pow<6>(_polar_y[_qp]) + Utility::pow<4>(_polar_y[_qp])*(-1 + Utility::pow<2>(_polar_z[_qp])) - Utility::pow<4>(_polar_z[_qp])*(1 + Utility::pow<2>(_polar_z[_qp])) + Utility::pow<4>(_polar_x[_qp])*(-1 + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])) + 
     Utility::pow<2>(_polar_y[_qp])*(Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp])) + Utility::pow<2>(_polar_x[_qp])*(Utility::pow<4>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]) + Utility::pow<4>(_polar_z[_qp]) + Utility::pow<2>(_polar_y[_qp])*(1 - 3*Utility::pow<2>(_polar_z[_qp])))));
}
