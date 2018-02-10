/**
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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "BulkEnergyDensity.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<BulkEnergyDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculates the free energy density dependent on the local polarization field.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1", "alpha1 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha11", "alpha11 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha12", "alpha12 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha111", "alpha111 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha112", "alpha112 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha123", "alpha123 coefficient of the Landau expansion");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

BulkEnergyDensity::BulkEnergyDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _alpha1(getParam<Real>("alpha1")),
  _alpha11(getParam<Real>("alpha11")),
  _alpha12(getParam<Real>("alpha12")),
  _alpha111(getParam<Real>("alpha111")),
  _alpha112(getParam<Real>("alpha112")),
  _alpha123(getParam<Real>("alpha123")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
BulkEnergyDensity::computeValue()
{
  return (
    _alpha1 * (Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) + std::pow(_polar_z[_qp] ,2))
  + _alpha11 * (Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_z[_qp]))
    + _alpha12 * (Utility::pow<2>(_polar_x[_qp]) * Utility::pow<2>(_polar_y[_qp])+
	      Utility::pow<2>(_polar_y[_qp]) * Utility::pow<2>(_polar_z[_qp])+
	      Utility::pow<2>(_polar_x[_qp]) * Utility::pow<2>(_polar_z[_qp]))+
    _alpha111 * (Utility::pow<6>(_polar_x[_qp]) + Utility::pow<6>(_polar_y[_qp]) + Utility::pow<6>(_polar_z[_qp]))+
    _alpha112 * (Utility::pow<4>(_polar_x[_qp]) * (Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))
	      + Utility::pow<4>(_polar_y[_qp]) * (Utility::pow<2>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp]))
	      + Utility::pow<4>(_polar_z[_qp]) * (Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp])))+
	  _alpha123 * (pow(_polar_x[_qp], 2) * Utility::pow<2>(_polar_y[_qp]) * Utility::pow<2>(_polar_z[_qp]))) * Utility::pow<3>(_len_scale);
}
