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

#include "TotalEnergySkFlow.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergySkFlow>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("Fbulk", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("Fwall", 0.0,  "name of wall energy postprocessor");
  params.addParam<PostprocessorName>("Felec", 0.0, "name of electrostatic energy postprocessor");
  params.addParam<PostprocessorName>("Fcoupled", 0.0, "name of the coupled energy postprocessor");
  params.addParam<PostprocessorName>("Fdepol", 0.0, "name of the depolarization energy postprocessor");
  return params;
}

TotalEnergySkFlow::TotalEnergySkFlow(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _Fbulk(getPostprocessorValue(getParam<PostprocessorName>("Fbulk"))),
  _Fwall(getPostprocessorValue(getParam<PostprocessorName>("Fwall"))),
  _Felec(getPostprocessorValue(getParam<PostprocessorName>("Felec"))),
  _Fcoupled(getPostprocessorValue(getParam<PostprocessorName>("Fcoupled"))),
  _Fdepol(getPostprocessorValue(getParam<PostprocessorName>("Fdepol")))
{
}

TotalEnergySkFlow::~TotalEnergySkFlow(){
}

void
TotalEnergySkFlow::initialize(){
}

void
TotalEnergySkFlow::execute(){
}

Real
TotalEnergySkFlow::getValue()
{
  return _Fbulk + _Fwall + _Felec + _Fcoupled + _Fdepol;
}
