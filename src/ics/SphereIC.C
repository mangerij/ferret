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

#include "SphereIC.h"
#include "Function.h"
#include "libmesh/point.h"

template<>
InputParameters validParams<SphereIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<FunctionName>("radial_function", "The function for radial");
  params.addRequiredParam<FunctionName>("polar_function", "The function for polar angle");
  params.addRequiredParam<FunctionName>("azimuthal_function", "The function for polar angle");
  params.addRequiredParam<unsigned int>("index", "which component to pick (0 for the first)");
  return params;
}



SphereIC::SphereIC(const InputParameters & parameters) :
  InitialCondition(parameters),
  _radial_func(getFunction("radial_function")),
  _polar_func(getFunction("polar_function")),
  _azimuthal_func(getFunction("azimuthal_function")),
  _index(getParam<unsigned int>("index"))
{
if(_index>2)
  mooseError("ERROR:index must be 0,1,2");
}

Real
SphereIC::value(const Point & p)
{
  Real radial=_radial_func.value(0.0,p);
  Real polar=_polar_func.value(0.0,p);
  Real azimuthal=_azimuthal_func.value(0.0,p);
  Real rv;
  switch(_index){
  case 0:
    rv=radial*cos(azimuthal)*cos(polar);
    break;
  case 1:
    rv=radial*cos(azimuthal)*sin(polar);
    break;
  case 2:
    rv=radial*sin(azimuthal);
    break;
  }
  return rv;
}
