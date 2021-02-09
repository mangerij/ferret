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

#include "BulkEnergyDensity.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", BulkEnergyDensity);

template<>
InputParameters validParams<BulkEnergyDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculates the free energy density dependent on the local polarization field.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  return params;
}

BulkEnergyDensity::BulkEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _alpha1(getMaterialProperty<Real>("alpha1")),
  _alpha11(getMaterialProperty<Real>("alpha11")),
  _alpha12(getMaterialProperty<Real>("alpha12")),
  _alpha111(getMaterialProperty<Real>("alpha111")),
  _alpha112(getMaterialProperty<Real>("alpha112")),
  _alpha123(getMaterialProperty<Real>("alpha123")),
  _alpha1111(getMaterialProperty<Real>("alpha1111")),
  _alpha1112(getMaterialProperty<Real>("alpha1112")),
  _alpha1122(getMaterialProperty<Real>("alpha1122")),
  _alpha1123(getMaterialProperty<Real>("alpha1123"))
{
}

Real
BulkEnergyDensity::computeValue()
{
  return (_alpha123[_qp]*Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + _alpha1[_qp]*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])) + _alpha12[_qp]*(Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])) + 
   _alpha11[_qp]*(Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_z[_qp])) + _alpha1123[_qp]*(Utility::pow<6>(_polar_x[_qp])*Utility::pow<2>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp])*Utility::pow<4>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*Utility::pow<4>(_polar_z[_qp])) + 
   _alpha1122[_qp]*(Utility::pow<4>(_polar_x[_qp])*Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_x[_qp])*Utility::pow<4>(_polar_z[_qp]) + Utility::pow<4>(_polar_y[_qp])*Utility::pow<4>(_polar_z[_qp])) + _alpha111[_qp]*(Utility::pow<6>(_polar_x[_qp]) + Utility::pow<6>(_polar_y[_qp]) + Utility::pow<6>(_polar_z[_qp])) + _alpha1111[_qp]*(Utility::pow<8>(_polar_x[_qp]) + Utility::pow<8>(_polar_y[_qp]) + Utility::pow<8>(_polar_z[_qp])) + 
   _alpha112[_qp]*((Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<4>(_polar_z[_qp]) + Utility::pow<4>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp])) + Utility::pow<4>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))) + 
   _alpha1112[_qp]*((Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<6>(_polar_z[_qp]) + Utility::pow<6>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp])) + Utility::pow<6>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))));
}
