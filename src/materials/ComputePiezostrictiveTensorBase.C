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

#include "ComputePiezostrictiveTensorBase.h"

InputParameters ComputePiezostrictiveTensorBase::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputePiezostrictiveTensorBase::ComputePiezostrictiveTensorBase(const InputParameters & parameters) :
    Material(parameters),
   _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
   _piezostrictive_tensor_name(_base_name + "piezostrictive_tensor"),
   _piezostrictive_tensor(declareProperty<RankThreeTensor>(_piezostrictive_tensor_name))
{
}

void
ComputePiezostrictiveTensorBase::computeQpProperties()
{
  computeQpPiezostrictiveTensor();
}
