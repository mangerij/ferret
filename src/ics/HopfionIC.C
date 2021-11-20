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

#include "HopfionIC.h"

#include "libmesh/point.h"
#include "libmesh/utility.h"
#include <cmath>

registerMooseObject("FerretApp", HopfionIC);

InputParameters HopfionIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addRequiredParam<Real>("Lhopf", "hopfion size value for the twirling pattern");
  params.addRequiredParam<Real>("p0", "spontaneous value");
  params.addRequiredParam<unsigned int>("component", "integer");
  return params;
}

HopfionIC::HopfionIC(const InputParameters & parameters) :
     InitialCondition(parameters),
    _Lhopf(getParam<Real>("Lhopf")),
    _p0(getParam<Real>("p0")),
    _component(getParam<unsigned int>("component"))
{
}

Real
HopfionIC::value(const Point & p)
{
  Real rho = std::sqrt(p(0)*p(0)+p(1)*p(1)+p(2)*p(2));
  Real Omg = std::tan(3.14159*p(2)/_Lhopf);
  Real a = (1.0/_Lhopf)*(1.0+(2*p(2)/_Lhopf)*(2*p(2)/_Lhopf))*(1.0/std::sin(3.14159*rho/(2*_Lhopf)));
  Real Lambda = (a*a*rho+Omg*Omg/4.0);
  Real th = std::atan2(p(1),p(0));

  if (_component == 0)
  {
    if (rho == 0.0)
    {
      return 0.0;
    }
    else
      return _p0*(4*rho*a*(-((-1 + Lambda)*std::cos(th)) + Omg*std::sin(th)))/Utility::pow<2>(1 + Lambda);
  }
  else if (_component == 1)
  {
    if (rho == 0.0)
    {
      return 0.0;
    }
    else
      return _p0*(4*rho*a*(-((-1 + Lambda)*std::sin(th)) + Omg*std::cos(th)))/Utility::pow<2>(1 + Lambda);
  }
  else if (_component == 2)
  {
    if (rho == 0.0)
    {
      return _p0;
    }
    else
      return _p0*(1.0 - 8.0*a*a*rho*rho/Utility::pow<2>(1 + Lambda));
  }
  else
    return 0.0;
}
