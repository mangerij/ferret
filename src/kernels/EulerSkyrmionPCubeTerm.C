/***************************************************************************/
/* This file is part of FERRET, an add-on module for MOOSE

/* FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

/* This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

/* You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

/****************************************************************************/

#include "EulerSkyrmionPCubeTerm.h"

class EulerSkyrmionPCubeTerm;

template<>
InputParameters validParams<EulerSkyrmionPCubeTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("P", "The polar magnitude variable");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("P0", 1.0, "the domain wall coefficient");
  return params;
}

EulerSkyrmionPCubeTerm::EulerSkyrmionPCubeTerm(const InputParameters & parameters)
  :Kernel(parameters),
   _P_var(coupled("P")),
   _P(coupledValue("P")),
   _len_scale(getParam<Real>("len_scale")),
   _P0(getParam<Real>("P0"))
{
}

Real
EulerSkyrmionPCubeTerm::computeQpResidual()
{
  return _test[_i][_qp] * _P[_qp] * _P[_qp] * _P[_qp] / (_P0 * _P0);
}

Real
EulerSkyrmionPCubeTerm::computeQpJacobian()
{
  return 2.0 * _test[_i][_qp] * _phi[_j][_qp] * _P[_qp] * _P[_qp] / (_P0 * _P0);
}
