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

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "EnergyRatePostprocessor.h"

registerMooseObject("FerretApp", EnergyRatePostprocessor);

template<>
InputParameters validParams<EnergyRatePostprocessor>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addClassDescription("Calculates the change of a postprocessor divided by the time step.");
  params.addRequiredParam<PostprocessorName>("postprocessor", "The name of the postprocessor used for exit criterion");
  params.addRequiredParam<PostprocessorName>("dt", "The dt postprocessor");
  return params;
}

EnergyRatePostprocessor::EnergyRatePostprocessor(const InputParameters & parameters) :
    GeneralPostprocessor(parameters),
    _postprocessor(getPostprocessorValue("postprocessor")),
    _postprocessor_old(getPostprocessorValueOld("postprocessor")),
    _dt(getPostprocessorValue("dt")),
    _dt_old(getPostprocessorValueOld("dt"))
{
}

void
EnergyRatePostprocessor::initialize(){
}

void
EnergyRatePostprocessor::execute(){
}

Real
EnergyRatePostprocessor::getValue()
{
  return fabs( ( fabs(_postprocessor)- fabs(_postprocessor_old) )*pow(fabs(_postprocessor),-1)) / _dt;
}
