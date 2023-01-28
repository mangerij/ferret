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

#include "InhomogeneousBulkEnergy.h"
#include "libmesh/utility.h"

registerMooseObject("FerretApp", InhomogeneousBulkEnergy);

InputParameters InhomogeneousBulkEnergy::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral whose integrand is the free energy density corresponding to the disordered materials coefficients.");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredCoupledVar("x", "The concentration");
  params.addParam<Real>("alpha01", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha011", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha0111", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha012", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha0112", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha0123", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha1111", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha1112", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha1122", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha1123", "The coefficients of the Landau expansion");
  params.addParam<Real>("b1", "The coupling coefficients to Sr concentration");
  params.addParam<Real>("b2", "The coupling coefficients to Sr concentration");
  params.addParam<Real>("b3", "The coupling coefficients to Sr concentration");
  params.addParam<Real>("T", "The temperature");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

InhomogeneousBulkEnergy::InhomogeneousBulkEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _x(coupledValue("x")),
   _alpha01(getParam<Real>("alpha01")),
   _alpha011(getParam<Real>("alpha011")),
   _alpha0111(getParam<Real>("alpha0111")),
   _alpha012(getParam<Real>("alpha012")),
   _alpha0112(getParam<Real>("alpha0112")),
   _alpha0123(getParam<Real>("alpha0123")),
   _alpha1111(getParam<Real>("alpha1111")),
   _alpha1112(getParam<Real>("alpha1112")),
   _alpha1122(getParam<Real>("alpha1122")),
   _alpha1123(getParam<Real>("alpha1123")),
   _b1(getParam<Real>("b1")),
   _b2(getParam<Real>("b2")),
   _b3(getParam<Real>("b3")),
   _T(getParam<Real>("T")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
InhomogeneousBulkEnergy::computeQpIntegral()
{
  return (Utility::pow<8>(_polar_x[_qp]) + Utility::pow<8>(_polar_y[_qp]) + Utility::pow<8>(_polar_z[_qp]))*_alpha1111 + ((Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]))*Utility::pow<6>(_polar_z[_qp]) + Utility::pow<6>(_polar_y[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_z[_qp])) + Utility::pow<6>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])))*_alpha1112 + 
   (Utility::pow<4>(_polar_y[_qp])*Utility::pow<4>(_polar_z[_qp]) + Utility::pow<4>(_polar_x[_qp])*(Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_z[_qp])))*_alpha1122 + Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])*(Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_alpha1123 + 
   (Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_z[_qp]))*_alpha011*(1 - _b2*_x[_qp]) + (Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp]) + Utility::pow<2>(_polar_x[_qp])*(Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp])))*_alpha012*(1 - _b2*_x[_qp]) + 
   ((Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_y[_qp]))*Utility::pow<2>(_polar_z[_qp]) + Utility::pow<2>(_polar_y[_qp])*(Utility::pow<4>(_polar_x[_qp]) + Utility::pow<4>(_polar_z[_qp])) + Utility::pow<2>(_polar_x[_qp])*(Utility::pow<4>(_polar_y[_qp]) + Utility::pow<4>(_polar_z[_qp])))*_alpha0112*(1 - _b3*_x[_qp]) + Utility::pow<2>(_polar_x[_qp])*Utility::pow<2>(_polar_y[_qp])*Utility::pow<2>(_polar_z[_qp])*_alpha0123*(1 - _b3*_x[_qp]) + 
   (Utility::pow<2>(_polar_x[_qp]) + Utility::pow<2>(_polar_y[_qp]) + Utility::pow<2>(_polar_z[_qp]))*_alpha01*(-2.5727418732060157 - _b1*_x[_qp] + ((std::exp(2.0*(160./_T))+1)/(std::exp(2.0*(160./_T))-1)));
}
