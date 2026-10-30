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

//* This file is adapted from part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LocalConservedNoiseBase.h"

#include "libmesh/quadrature.h"

InputParameters
LocalConservedNoiseBase::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
  params.addParam<bool>("conserved", true, "Subtract the volume average so the noise integrates to zero (concentration-type fields); false keeps the uniform mode (non-conserved order parameters, macrospins)");
  params.addParam<bool>("white_noise", false, "Divide each value by sqrt(JxW dt): space-time white noise with unit strength, so a source -test*sqrt(2D)*value has the fluctuation-dissipation covariance 2D M/dt (M the consistent mass matrix)");
  return params;
}

LocalConservedNoiseBase::LocalConservedNoiseBase(const InputParameters & parameters)
  : LocalConservedNoiseInterface(parameters),
    _conserved(getParam<bool>("conserved")),
    _white_noise(getParam<bool>("white_noise"))
{
}

void
LocalConservedNoiseBase::initialize()
{
  _random_data.clear();
  _integral = 0.0;
  _volume = 0.0;
}

void
LocalConservedNoiseBase::execute()
{
  // reserve space for each quadrature point in the element
  std::vector<Real> & me = _random_data[_current_elem->id()] =
      std::vector<Real>(_qrule->n_points());

  // store a random number for each quadrature point
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    me[_qp] = getQpRandom();
    if (_white_noise)
      me[_qp] = _dt > 0.0 ? me[_qp] / std::sqrt(_JxW[_qp] * _coord[_qp] * _dt) : 0.0;
    _integral += _JxW[_qp] * _coord[_qp] * me[_qp];
    _volume += _JxW[_qp] * _coord[_qp];
  }
}

void
LocalConservedNoiseBase::threadJoin(const UserObject & y)
{
  const LocalConservedNoiseBase & uo = static_cast<const LocalConservedNoiseBase &>(y);

  _random_data.insert(uo._random_data.begin(), uo._random_data.end());
  _integral += uo._integral;
  _volume += uo._volume;
}

void
LocalConservedNoiseBase::finalize()
{
  gatherSum(_integral);
  gatherSum(_volume);

  _offset = _conserved ? _integral / _volume : 0.0;
}

Real
LocalConservedNoiseBase::getQpValue(dof_id_type element_id, unsigned int qp) const
{
  const auto it_pair = _random_data.find(element_id);

  if (it_pair == _random_data.end())
    mooseError("Element not found.");
  else
  {
    libmesh_assert_less(qp, it_pair->second.size());
    return it_pair->second[qp] - _offset;
  }
}
