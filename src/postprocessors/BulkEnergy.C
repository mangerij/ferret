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

#include "BulkEnergy.h"

constexpr Real BulkEnergy::_default_uniform_val;

template<>
InputParameters validParams<BulkEnergy>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha3", BulkEnergy::_default_uniform_val, "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha11", "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha33", BulkEnergy::_default_uniform_val, "The coefficients of the Landau expansion");
  params.addParam<Real>("alpha13", BulkEnergy::_default_uniform_val, "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha12", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha111", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha112", "The coefficients of the Landau expansion");
  params.addRequiredParam<Real>("alpha123", "The coefficients of the Landau expansion");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

BulkEnergy::BulkEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _alpha1(getParam<Real>("alpha1")),
  _alpha3(getParam<Real>("alpha3") == _default_uniform_val ? _alpha1 : getParam<Real>("alpha3")),
  _alpha11(getParam<Real>("alpha11")),
  _alpha33(getParam<Real>("alpha33") == _default_uniform_val ? _alpha11 : getParam<Real>("alpha33")),
  _alpha12(getParam<Real>("alpha12")),
  _alpha13(getParam<Real>("alpha13") == _default_uniform_val ? _alpha12 : getParam<Real>("alpha13")),
  _alpha111(getParam<Real>("alpha111")),
  _alpha112(getParam<Real>("alpha112")),
  _alpha123(getParam<Real>("alpha123")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
BulkEnergy::computeQpIntegral()
{
  return (
        _alpha1 * ( std::pow(_polar_x[_qp], 2.0) + std::pow(_polar_y[_qp], 2.0) ) + 
        _alpha3 * std::pow(_polar_z[_qp], 2.0) + 
        _alpha11 * ( std::pow(_polar_x[_qp], 4.0) + std::pow(_polar_y[_qp], 4.0) ) +
        _alpha33 * std::pow(_polar_z[_qp], 4.0) + 
        _alpha13 * ( std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) + std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0) ) +
        _alpha12 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0) +
        _alpha111 * (std::pow(_polar_x[_qp], 6.0) + std::pow(_polar_y[_qp], 6.0) + std::pow(_polar_z[_qp], 6.0)) +
        _alpha112 * (
            std::pow(_polar_x[_qp], 4.0) * (std::pow(_polar_y[_qp], 2.0) + std::pow(_polar_z[_qp], 2.0)) + 
            std::pow(_polar_y[_qp], 4.0) * (std::pow(_polar_x[_qp], 2.0) + std::pow(_polar_z[_qp], 2.0)) + 
            std::pow(_polar_z[_qp], 4.0) * (std::pow(_polar_y[_qp], 2.0) + std::pow(_polar_x[_qp], 2.0)) 
        ) + 
        _alpha123 * std::pow(_polar_x[_qp], 2.0) * std::pow(_polar_y[_qp], 2.0) * std::pow(_polar_z[_qp], 2.0)
    );
}
