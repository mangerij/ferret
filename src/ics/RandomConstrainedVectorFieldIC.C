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

#include "RandomConstrainedVectorFieldIC.h"

#include "libmesh/point.h"
#include <cmath>

registerMooseObject("FerretApp", RandomConstrainedVectorFieldIC);

InputParameters RandomConstrainedVectorFieldIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this IC is set. (0 for x, 1 for y, 2.0 for z)");
  params.addRequiredCoupledVar("phi", "The field of random azimuthal angles");
  params.addRequiredCoupledVar("theta", "The field of random polar angles");
  params.addRequiredParam<Real>("M0s", "M0s");
  return params;
}

RandomConstrainedVectorFieldIC::RandomConstrainedVectorFieldIC(const InputParameters & parameters) :
    InitialCondition(parameters),
    _component(getParam<unsigned int>("component")),
    _phi(coupledValue("phi")),
    _theta(coupledValue("theta")),
    _M0s(getParam<Real>("M0s"))
{
}

Real
RandomConstrainedVectorFieldIC::value(const Point & /* p */)
{
  if (_component == 0)
  {
    return _M0s*std::cos(_phi[_qp]) * std::sin(_theta[_qp]);
  }
  else if (_component == 1)
  {
    return _M0s*std::sin(_phi[_qp]) * std::sin(_theta[_qp]);
  }
  else if (_component == 2)
  {
    return _M0s*std::cos(_theta[_qp]);
  }
  else
    return 0.0;
}
