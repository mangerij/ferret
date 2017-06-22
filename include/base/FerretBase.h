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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

#ifndef FERRETBASE_H
#define FERRETBASE_H

#include <string>
#include <set>
#include "InputParameters.h"
#include "MooseObject.h"

class FerretBase
{
 public:
  FerretBase(InputParameters parameters);
 protected:
  std::vector<std::string> _debug_vec;
  std::string              _class;
  std::set<std::string>    _debug_set;
  bool debug(const std::string&);
};

#endif //FERRETBASE_H
