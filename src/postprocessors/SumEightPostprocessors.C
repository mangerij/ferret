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

#include "SumEightPostprocessors.h"
#include<iostream>
template<>
InputParameters validParams<SumEightPostprocessors>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<PostprocessorName>("var0", "name of first postprocessor");
  params.addRequiredParam<PostprocessorName>("var1", "name of second postprocessor");
  params.addRequiredParam<PostprocessorName>("var2", "name of third postprocessor");
  params.addRequiredParam<PostprocessorName>("var3", "name of fourth postprocessor");
  params.addRequiredParam<PostprocessorName>("var4", "name of fifth postprocessor");
  params.addRequiredParam<PostprocessorName>("var5", "name of sixth postprocessor");
  params.addRequiredParam<PostprocessorName>("var6", "name of seventh postprocessor");
  params.addRequiredParam<PostprocessorName>("var7", "name of eighth postprocessor");
  return params;
}

SumEightPostprocessors::SumEightPostprocessors(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
  _var0(getPostprocessorValue("var0")),
  _var1(getPostprocessorValue("var1")),
  _var2(getPostprocessorValue("var2")),
  _var3(getPostprocessorValue("var3")),
  _var4(getPostprocessorValue("var4")),
  _var5(getPostprocessorValue("var5")),
  _var6(getPostprocessorValue("var6")),
  _var7(getPostprocessorValue("var7"))
{
}

SumEightPostprocessors::~SumEightPostprocessors(){
}

void
SumEightPostprocessors::initialize(){
}

void
SumEightPostprocessors::execute(){
}

Real
SumEightPostprocessors::getValue()
{
  return _var0 + _var1 + _var2 + _var3 + _var4 + _var5 + _var6 + _var7;
}
