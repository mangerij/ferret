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

#include "FerretBase.h"

template<>
InputParameters validParams<FerretBase>()
{
  InputParameters params = emptyInputParameters();
  params.addParam<std::vector<std::string> >("Debug", std::vector<std::string>(), "Debug strings");
  params.addParam<std::string>("class","","Class name");
  return params;
}


FerretBase::FerretBase(InputParameters parameters):
  _debug_vec(parameters.get<std::vector<std::string> >("Debug")),
  _class(parameters.get<std::string>("class"))
{
  /// std::copy(_debug_vec.begin(),_debug_vec.end(),std::insert_iterator<std::set<std::string> >(_debug_set,_debug_set.begin()));
  for (unsigned int i = 0; i < _debug_vec.size(); ++i) {
    _debug_set.insert(_debug_vec[i]);
  }
}

bool
FerretBase::debug(const std::string& str)
{
  std::string s;
  if (_class.size()) s = _class+"::"+str; else s = str;
  unsigned int count = _debug_set.count(s);
  bool res = (count > 0);
  return res;
}
