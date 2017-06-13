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

#include "AcceptorIonContribution.h"

class AcceptorIonContribution;

template<>
InputParameters validParams<AcceptorIonContribution>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<Real>("q", 1.0, "The charge");
  params.addParam<Real>("Na", 1.0, "The ionized acceptors of the problem");
  params.addParam<Real>("len_scale", 1.0, "The length scale of the unit");
  return params;
}

AcceptorIonContribution::AcceptorIonContribution(const InputParameters & parameters)
  :Kernel(parameters),
   _q(getParam<Real>("q")),
   _Na(getParam<Real>("Na")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
AcceptorIonContribution::computeQpResidual()
{
  return std::pow(_len_scale, 2.0) * _q * _Na * _test[_i][_qp];
}

Real
AcceptorIonContribution::computeQpJacobian()
{
  return 0.0;
}
//need off diagon for c
