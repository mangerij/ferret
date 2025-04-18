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

#include "DepolarizationEnergy.h"

registerMooseObject("FerretApp", DepolarizationEnergy);

InputParameters DepolarizationEnergy::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over a fictious depolarization field energy density");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("avePz", 1.0, "the average polarization due to the depol field term");
  params.addParam<Real>("lambda", 1.0, "the screening length term");
  params.addParam<Real>("permitivitty", 1.0, "the permitivitty term");
  return params;
}

DepolarizationEnergy::DepolarizationEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _polar_z(coupledValue("polar_z")),
   _len_scale(getParam<Real>("len_scale")),
   _avePz(getParam<Real>("avePz")),
   _lambda(getParam<Real>("lambda")),
   _permitivitty(getParam<Real>("permitivitty"))
{
}

Real
DepolarizationEnergy::computeQpIntegral()
{
  return 0.5 * _lambda * (1.0 / _permitivitty) * _polar_z[_qp] * _avePz * std::pow(_len_scale, 3.0);
}
