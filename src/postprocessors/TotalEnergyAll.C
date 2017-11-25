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

#include "TotalEnergyAll.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergyAll>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("Fbulk", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("Fwall", 0.0,  "name of wall energy postprocessor");
  params.addParam<PostprocessorName>("Felec", 0.0, "name of electrostatic energy postprocessor");
  params.addParam<PostprocessorName>("Fcoupled", 0.0, "name of the coupled energy postprocessor");
  params.addParam<PostprocessorName>("Felastic", 0.0, "name of the elastic energy postprocessor");
  params.addParam<PostprocessorName>("Fextelec", 0.0, "name of the external electrostatic energy postprocessor");
  return params;
}

TotalEnergyAll::TotalEnergyAll(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _Fbulk(getPostprocessorValue(getParam<PostprocessorName>("Fbulk"))),
  _Fwall(getPostprocessorValue(getParam<PostprocessorName>("Fwall"))),
  _Felec(getPostprocessorValue(getParam<PostprocessorName>("Felec"))),
  _Fcoupled(getPostprocessorValue(getParam<PostprocessorName>("Fcoupled"))),
  _Felastic(getPostprocessorValue(getParam<PostprocessorName>("Felastic"))),
  _Fextelec(getPostprocessorValue(getParam<PostprocessorName>("Fextelec")))
{
}

TotalEnergyAll::~TotalEnergyAll(){
}

void
TotalEnergyAll::initialize(){
}

void
TotalEnergyAll::execute(){
}

Real
TotalEnergyAll::getValue()
{
  return _Fbulk + _Fwall + _Felec + _Fcoupled + _Felastic + _Fextelec;
}
