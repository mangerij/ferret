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

#include "PercentChangePostprocessor.h"

template<>
InputParameters validParams<PercentChangePostprocessor>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addClassDescription("Calculates the absolute value change between adjacent postprocessors in time");
  params.addRequiredParam<PostprocessorName>("postprocessor", "The name of the postprocessor used for exit criterion");
  return params;
}

PercentChangePostprocessor::PercentChangePostprocessor(const InputParameters & parameters) :
    GeneralPostprocessor(parameters),
    _postprocessor(getPostprocessorValue("postprocessor")),
    _postprocessor_old(getPostprocessorValueOld("postprocessor"))

{
}

void
PercentChangePostprocessor::initialize(){
}

void
PercentChangePostprocessor::execute(){
}

Real
PercentChangePostprocessor::getValue()
{
  return fabs( ( fabs(_postprocessor)- fabs(_postprocessor_old) )*pow(fabs(_postprocessor),-1));
}
