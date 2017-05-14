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

#include "EulerSkyrmionPDepolTerm.h"

class EulerSkyrmionPDepolTerm;

template<>
InputParameters validParams<EulerSkyrmionPDepolTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("theta", "The theta variable");
  params.addParam<Real>("edep", 1.0, "the depolarization field strength");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

EulerSkyrmionPDepolTerm::EulerSkyrmionPDepolTerm(const InputParameters & parameters)
  :Kernel(parameters),
   _theta_var(coupled("theta")),
   _theta(coupledValue("theta")),
   _edep(getParam<Real>("edep")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
EulerSkyrmionPDepolTerm::computeQpResidual()
{
  return - _test[_i][_qp] * _edep * std::cos(_theta[_qp]);
}

Real
EulerSkyrmionPDepolTerm::computeQpJacobian()
{
  return 0.0;
}

Real
EulerSkyrmionPDepolTerm::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _theta_var)
    return _test[_i][_qp] * _phi[_j][_qp] * _edep * std::sin(_theta[_qp]);
  else
    return 0.0;
}
