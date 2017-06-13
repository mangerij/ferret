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

#include "TotalEnergyPSTOcoupled.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergyPSTOcoupled>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("FbulkPSTO", 0.0, "name of bulk energy postprocessor");
  params.addParam<PostprocessorName>("Fwall", 0.0,  "name of wall energy postprocessor");
  params.addParam<PostprocessorName>("Felec", 0.0,  "name of electrical energy postprocessor");
  params.addParam<PostprocessorName>("Fcoupled", 0.0,  "name of coupled energy postprocessor");

  return params;
}

TotalEnergyPSTOcoupled::TotalEnergyPSTOcoupled(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _FbulkPSTO(getPostprocessorValue(getParam<PostprocessorName>("FbulkPSTO"))),
  _Fwall(getPostprocessorValue(getParam<PostprocessorName>("Fwall"))),
  _Felec(getPostprocessorValue(getParam<PostprocessorName>("Felec"))),
  _Fcoupled(getPostprocessorValue(getParam<PostprocessorName>("Fcoupled")))


{
}

TotalEnergyPSTOcoupled::~TotalEnergyPSTOcoupled(){
}

void
TotalEnergyPSTOcoupled::initialize(){
}

void
TotalEnergyPSTOcoupled::execute(){
}

Real
TotalEnergyPSTOcoupled::getValue()
{
  ///  return _bulk_energy + _wall_energy + _electrical_energy + _coupled_energy +;
  return _FbulkPSTO + _Fwall + _Felec + _Fcoupled;
}
