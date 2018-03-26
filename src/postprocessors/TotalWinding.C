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

#include "TotalWinding.h"

registerMooseObject("FerretApp", TotalWinding);

template<>
InputParameters validParams<TotalWinding>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculates an integral over winding number density.");
  params.addRequiredCoupledVar("q", "The winding number density");
  return params;
}

TotalWinding::TotalWinding(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _q(coupledValue("q"))
{
}

Real
TotalWinding::computeQpIntegral()
{
  return _q[_qp];
}
