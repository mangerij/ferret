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

#include "GrainSize.h"
#include<iostream>

template<>
InputParameters validParams<GrainSize>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("vol", 0.0, "name of volume postprocessor");
  return params;
}

GrainSize::GrainSize(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _vol(getPostprocessorValue(getParam<PostprocessorName>("vol")))
{
}

GrainSize::~GrainSize(){
}

void
GrainSize::initialize(){
}

void
GrainSize::execute(){
}

Real
GrainSize::getValue()
{
  return 2.0 * std::pow(3.0 * _vol / (4.0 * 3.14159), 0.33333);
}
