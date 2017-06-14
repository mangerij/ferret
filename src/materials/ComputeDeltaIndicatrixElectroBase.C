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

#include "ComputeDeltaIndicatrixElectroBase.h"

template<>
InputParameters validParams<ComputeDeltaIndicatrixElectroBase>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}


ComputeDeltaIndicatrixElectroBase::ComputeDeltaIndicatrixElectroBase(const InputParameters & parameters) :
    Material(parameters),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
   _delta_indicatrix_electro_name(_base_name + "delta_indicatrix"),
   _delta_indicatrix_electro(declareProperty<RankTwoTensor>("delta_indicatrix_electro"))
{
}

void
ComputeDeltaIndicatrixElectroBase::computeQpProperties()
{
  computeQpDeltaIndicatrixElectro();
}
