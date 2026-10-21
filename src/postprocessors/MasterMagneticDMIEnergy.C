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

#include "MasterMagneticDMIEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", MasterMagneticDMIEnergy);

InputParameters MasterMagneticDMIEnergy::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over the bulk or interfacial Dzyaloshinskii-Moriya energy density of a ferromagnet.");
  params.addRequiredCoupledVar("mag_x", "The x component of the magnetization");
  params.addRequiredCoupledVar("mag_y", "The y component of the magnetization");
  params.addCoupledVar("mag_z", 0.0, "The z component of the magnetization");
  params.addParam<MooseEnum>("dmi_type", MooseEnum("bulk interfacial", "bulk"), "bulk: w = D m.(curl m); interfacial: w = D [mz div(m) - (m.grad) mz], interface normal +z");
  params.addParam<Real>("energy_scale", 1.0, "the energy scale, useful for transition between eV and aJ");
  return params;
}

MasterMagneticDMIEnergy::MasterMagneticDMIEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _bulk(getParam<MooseEnum>("dmi_type") == "bulk"),
  _mag_x(coupledValue("mag_x")),
  _mag_y(coupledValue("mag_y")),
  _mag_z(coupledValue("mag_z")),
  _mag_x_grad(coupledGradient("mag_x")),
  _mag_y_grad(coupledGradient("mag_y")),
  _mag_z_grad(coupledGradient("mag_z")),
  _D(getMaterialProperty<Real>("D")),
  _energy_scale(getParam<Real>("energy_scale"))
{
}

Real
MasterMagneticDMIEnergy::computeQpIntegral()
{
  if (_bulk)
    return _energy_scale*_D[_qp]*(-_mag_x_grad[_qp](1)*_mag_z[_qp] + _mag_x_grad[_qp](2)*_mag_y[_qp] + _mag_y_grad[_qp](0)*_mag_z[_qp] - _mag_y_grad[_qp](2)*_mag_x[_qp] - _mag_z_grad[_qp](0)*_mag_y[_qp] + _mag_z_grad[_qp](1)*_mag_x[_qp]);
  else
    return _energy_scale*_D[_qp]*(_mag_x_grad[_qp](0)*_mag_z[_qp] + _mag_y_grad[_qp](1)*_mag_z[_qp] - _mag_z_grad[_qp](0)*_mag_x[_qp] - _mag_z_grad[_qp](1)*_mag_y[_qp]);
}
