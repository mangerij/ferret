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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#include "FluctuationsIC.h"

#include "libmesh/point.h"
#include <cmath>

template<>
InputParameters validParams<FluctuationsIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("epsilon", "epsilon value for psudorandom function");
  params.addParam<Real>("base_value", 0.0, "Value that is fluctuated about");
  params.addRequiredParam<Point>("q1", "first q-vector for pseudorandom function");
  params.addRequiredParam<Point>("q2", "second q-vector for pseudorandom function");
  params.addRequiredParam<Point>("q3", "third q-vector for pseudorandom function");
  params.addRequiredParam<Point>("q4", "fourth q-vector for pseudorandom function");
  params.addRequiredParam<Real>("h", "h factor on cos^2 term");

  return params;
}

FluctuationsIC::FluctuationsIC(const InputParameters & parameters) :
    InitialCondition(parameters),
    _epsilon(getParam<Real>("epsilon")),
    _base_value(getParam<Real>("base_value")),
    _q1(getParam<Point>("q1")),
    _q2(getParam<Point>("q2")),
    _q3(getParam<Point>("q3")),
    _q4(getParam<Point>("q4")),
    _h(getParam<Real>("h"))
{
}

Real
FluctuationsIC::value(const Point & p)
{
  Real delta_eta = (std::cos(_q1(0) * p(0) - 4) * std::cos(_q1(1) * p(1))
    + std::cos(_q2(0)*p(0)) * std::cos(_q2(1) * p(1))
    + _h * std::pow( std::cos(_q1(0) * p(1) + _q3(1) * p(1)) * std::cos(_q4(0) * p(0) - _q4(1) * p(1)) , 2.0 ) - std::cos(_q2(0) * p(0) - 4) * std::cos(_q1(1) * p(1))
    + std::cos(_q2(1)*p(2)) * std::cos(_q2(1) * p(0))
    + _h * std::pow( std::cos(_q3(0) * p(0) + _q3(0) * p(2)) * std::cos(_q4(1) * p(1) - _q4(0) * p(2)) , 2.0 ) );

  return _epsilon * std::pow( delta_eta, 2.0 ) - _epsilon;
}
