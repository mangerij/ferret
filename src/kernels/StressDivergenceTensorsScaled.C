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

#include "StressDivergenceTensorsScaled.h"

template<>
InputParameters validParams<StressDivergenceTensorsScaled>()
{
  InputParameters params = validParams<StressDivergenceTensors>();
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

StressDivergenceTensorsScaled::StressDivergenceTensorsScaled(const InputParameters & parameters)
  :StressDivergenceTensors(parameters),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
StressDivergenceTensorsScaled::computeQpResidual()
{
  return std::pow(_len_scale, 2.0) * StressDivergenceTensors::computeQpResidual();
}

Real
StressDivergenceTensorsScaled::computeQpJacobian()
{
  return  std::pow(_len_scale, 2.0) * StressDivergenceTensors::computeQpJacobian();
}

Real
StressDivergenceTensorsScaled::computeQpOffDiagJacobian(unsigned int jvar)
{
  return  std::pow(_len_scale, 2.0) * StressDivergenceTensors::computeQpOffDiagJacobian(jvar);
}
