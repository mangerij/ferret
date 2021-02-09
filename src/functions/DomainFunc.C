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

#include "DomainFunc.h"
#include "MooseRandom.h"

registerMooseObject("FerretApp", DomainFunc);

template <>
InputParameters
validParams<DomainFunc>()
{
  InputParameters params = validParams<Function>();
  params.addRequiredParam<Real>("ax", "Width of ferroelectric");
  params.addRequiredParam<Real>("af", "Thickness of ferroelectric");
  params.addRequiredParam<Real>("min", "Min value for random perturbation of polarization");
  params.addRequiredParam<Real>("max", "Max value for random perturbation of polarization");
  return params;
}

DomainFunc::DomainFunc(const InputParameters & parameters)
  : Function(parameters),
  _ax(getParam<Real>("ax")),
  _af(getParam<Real>("af")),
  _min(getParam<Real>("min")),
  _max(getParam<Real>("max"))
{
    MooseRandom::seed(0.0);
}

Real
DomainFunc::value(Real /*t*/, const Point & p) const
{
  Real rand_num = MooseRandom::rand();

  // Between 0 and range
  rand_num *= (_max - _min);

  // Between min and max
  rand_num += _min;

  return 1.4 * std::cos(2.0 * libMesh::pi * p(0) / _ax) * std::cos(libMesh::pi * p(2) / _af) + rand_num;
}
