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

#include "AdhocConstIC.h"

#include "libmesh/point.h"
#include <limits>
#include <cmath>
template<>
InputParameters validParams<AdhocConstIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("value0", "value for z<0.5");
  params.addRequiredParam<Real>("value1", "value for z>=0.5");
  return params;
}



AdhocConstIC::AdhocConstIC(const InputParameters & parameters) :
  InitialCondition(parameters),
  _val0(getParam<Real>("value0")),
  _val1(getParam<Real>("value1"))
{
}

Real
AdhocConstIC::value(const Point & p)
{
  if(p(2)-0.5<1000*std::numeric_limits<Real>::epsilon()) return _val0;
  else return _val1;
}
