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

#include "TotalEnergyBFOElec.h"
#include<iostream>
template<>
InputParameters validParams<TotalEnergyBFOElec>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addParam<PostprocessorName>("FbP", 0.0, "name of polar bulk energy postprocessor");
  params.addParam<PostprocessorName>("FbA", 0.0,  "name of roto bulk energy postprocessor");
  params.addParam<PostprocessorName>("FgP", 0.0,  "name of polar gradient energy postprocessor");
  params.addParam<PostprocessorName>("FgA", 0.0,  "name of roto gradient energy postprocessor");
  params.addParam<PostprocessorName>("FcPA", 0.0,  "name of coupled rotopolar energy postprocessor");
  params.addParam<PostprocessorName>("FcPu", 0.0,  "name of coupled electrostrictive energy postprocessor");
  params.addParam<PostprocessorName>("FcAu", 0.0,  "name of coupled rotostrictive energy postprocessor");
  params.addParam<PostprocessorName>("Felu", 0.0,  "name of elastic energy postprocessor");
  params.addParam<PostprocessorName>("Felec", 0.0,  "name of electrostatic energy postprocessor");
  return params;
}

TotalEnergyBFOElec::TotalEnergyBFOElec(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _FbP(getPostprocessorValue(getParam<PostprocessorName>("FbP"))),
  _FbA(getPostprocessorValue(getParam<PostprocessorName>("FbA"))),
  _FgP(getPostprocessorValue(getParam<PostprocessorName>("FgP"))),
  _FgA(getPostprocessorValue(getParam<PostprocessorName>("FgA"))),
  _FcPA(getPostprocessorValue(getParam<PostprocessorName>("FcPA"))),
  _FcPu(getPostprocessorValue(getParam<PostprocessorName>("FcPu"))),
  _FcAu(getPostprocessorValue(getParam<PostprocessorName>("FcAu"))),
  _Felu(getPostprocessorValue(getParam<PostprocessorName>("Felu"))),
  _Felec(getPostprocessorValue(getParam<PostprocessorName>("Felec")))
{
}

TotalEnergyBFOElec::~TotalEnergyBFOElec(){
}

void
TotalEnergyBFOElec::initialize(){
}

void
TotalEnergyBFOElec::execute(){
}

Real
TotalEnergyBFOElec::getValue()
{
  return _FbP + _FbA + _FgP + _FgA + _FcPA + _FcPu + _FcAu + _Felu + _Felec;
}
