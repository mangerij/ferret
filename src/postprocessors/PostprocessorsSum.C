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

#include "PostprocessorsSum.h"
#include<iostream>
template<>
InputParameters validParams<PostprocessorsSum>()
{

  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<unsigned int>("num_pp", "How many postprocessors to sum? Need to be 2 < num_pp < 9.");
  params.addRequiredParam<PostprocessorName>("var0", "name of first postprocessor");
  params.addRequiredParam<PostprocessorName>("var1", "name of second postprocessor");
  params.addParam<PostprocessorName>("var2", 0.0, "name of third postprocessor");
  params.addParam<PostprocessorName>("var3", 0.0, "name of fourth postprocessor");
  params.addParam<PostprocessorName>("var4", 0.0, "name of fifth postprocessor");
  params.addParam<PostprocessorName>("var5", 0.0, "name of sixth postprocessor");
  params.addParam<PostprocessorName>("var6", 0.0, "name of seventh postprocessor");
  params.addParam<PostprocessorName>("var7", 0.0, "name of eighth postprocessor");
  params.addParam<PostprocessorName>("var8", 0.0, "name of ninth postprocessor");
  return params;
}

PostprocessorsSum::PostprocessorsSum(const InputParameters & parameters) :
  GeneralPostprocessor(parameters),
   _num_pp(getParam<unsigned int>("num_pp")),
  _var0(getPostprocessorValue("var0")),
  _var1(getPostprocessorValue("var1")),
  _var2(getPostprocessorValue("var2")),
  _var3(getPostprocessorValue("var3")),
  _var4(getPostprocessorValue("var4")),
  _var5(getPostprocessorValue("var5")),
  _var6(getPostprocessorValue("var6")),
  _var7(getPostprocessorValue("var7")),
  _var8(getPostprocessorValue("var8"))
{
}

PostprocessorsSum::~PostprocessorsSum(){
}

void
PostprocessorsSum::initialize(){
}

void
PostprocessorsSum::execute(){
}

Real
PostprocessorsSum::getValue()
{
  if (_num_pp == 2)
  {
    return _var0 + _var1;
  }
  else if (_num_pp == 3)
  {
    return _var0 + _var1 + _var2;
  }
  else if (_num_pp == 4)
  {
    return _var0 + _var1 + _var2 + _var3;
  }
  else if (_num_pp == 5)
  {
    return _var0 + _var1 + _var2 + _var3 + _var4;
  }
  else if (_num_pp == 6)
  {
    return _var0 + _var1 + _var2 + _var3 + _var4 + _var5;
  }
  else if (_num_pp == 7)
  {
    return _var0 + _var1 + _var2 + _var3 + _var4 + _var5 + _var6;
  }
  else if (_num_pp == 8)
  {
    return _var0 + _var1 + _var2 + _var3 + _var4 + _var5 + _var6 + _var7;
  }
  else if (_num_pp == 9)
  {
    return _var0 + _var1 + _var2 + _var3 + _var4 + _var5 + _var6 + _var7 + _var8;
  }
  else
    return 0.0;
}
