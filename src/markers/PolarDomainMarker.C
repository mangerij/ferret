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

#include "PolarDomainMarker.h"
#include "FEProblem.h"
#include "MooseEnum.h"
#include "RankTwoTensor.h"
#include <math.h>

registerMooseObject("FerretApp", PolarDomainMarker);

InputParameters
PolarDomainMarker::validParams()
{
  InputParameters params = QuadraturePointMarker::validParams();
  params.addRequiredParam<Real>("lower_bound", "The lower bound value for the range.");
  params.addRequiredParam<Real>("upper_bound", "The upper bound value for the range.");
  params.addParam<Real>("buffer_size",
                        0.0,
                        "A buffer zone value added to both ends of the range "
                        "where a third_state marker can be returned.");
  params.addClassDescription("Mark elements for adaptivity based on the supplied upper and lower "
                             "bounds and the specified variable.");
  return params;
}

PolarDomainMarker::PolarDomainMarker(const InputParameters & parameters)
  : QuadraturePointMarker(parameters),
    _lower_bound(parameters.get<Real>("lower_bound")),
    _upper_bound(parameters.get<Real>("upper_bound")),
    _buffer_size(parameters.get<Real>("buffer_size")),
    _inside(getParam<bool>("invert") ? COARSEN : DO_NOTHING),
    _outside(getParam<bool>("invert") ? DO_NOTHING : COARSEN)
{
  if (_upper_bound < _lower_bound)
    mooseError("Invalid bounds specified (upper_bound < lower_bound)");

  if (_buffer_size < 0.0)
    mooseError("Buffer size must be non-negative: ", _buffer_size);
}

Marker::MarkerValue
PolarDomainMarker::computeQpMarker()
{
  // Is the variable value inside the range?
  if (abs(_u[_qp]) >= _lower_bound && abs(_u[_qp]) <= _upper_bound)
    return _inside;

  // How about the buffer zone?
  if (abs(_u[_qp]) >= _lower_bound - _buffer_size && abs(_u[_qp]) <= _upper_bound + _buffer_size)
    return _third_state;

  // Must be outside the range
  return _outside;
}
