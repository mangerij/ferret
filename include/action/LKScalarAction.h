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

   For help with FERRET please contact J. Mangeri <johnma@dtu.dk>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#pragma once

#include "Action.h"
#include "libmesh/fe_type.h"

class LKScalarAction : public Action
{
public:
  LKScalarAction(const InputParameters & params);
  static InputParameters validParams();
  virtual void act() override;

protected:
  /// Equation type, transient or steady-state
  //MooseEnum _type;

  /// FE type for various variables
  FEType _fe_type;

};
