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

//* This file is adapted from part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "MultiAppAddTransfer.h"
#include "FEProblemBase.h"
#include "MultiApp.h"

registerMooseObject("MooseApp", MultiAppAddTransfer);

InputParameters
MultiAppAddTransfer::validParams()
{
  InputParameters params = MultiAppFieldAddTransfer::validParams();
  params.addRequiredParam<std::vector<AuxVariableName>>(
      "variable", "The auxiliary variable to store the transferred values in.");
  params.addRequiredParam<std::vector<VariableName>>("source_variable",
                                                     "The variable to transfer from.");

  params.addClassDescription(
      "Adds variables (nonlinear and auxiliary) between multiapps that have identical meshes.");
  return params;
}

MultiAppAddTransfer::MultiAppAddTransfer(const InputParameters & parameters)
  : MultiAppFieldAddTransfer(parameters),
    _from_var_names(getParam<std::vector<VariableName>>("source_variable")),
    _to_var_names(getParam<std::vector<AuxVariableName>>("variable"))
{
  /* Right now, most of transfers support one variable only */
  _to_var_name = _to_var_names[0];
  _from_var_name = _from_var_names[0];
}

void
MultiAppAddTransfer::execute()
{
  _console << "Beginning MultiAppAddTransfer " << name() << std::endl;

  if (_current_direction == TO_MULTIAPP)
  {
    FEProblemBase & from_problem = getToMultiApp()->problemBase();
    for (unsigned int i = 0; i < getToMultiApp()->numGlobalApps(); i++)
      if (getToMultiApp()->hasLocalApp(i))
        transfer(getToMultiApp()->appProblemBase(i), from_problem);
  }

  else if (_current_direction == FROM_MULTIAPP)
  {
    FEProblemBase & to_problem = getFromMultiApp()->problemBase();
    for (unsigned int i = 0; i < getFromMultiApp()->numGlobalApps(); i++)
      if (getFromMultiApp()->hasLocalApp(i))
        transfer(to_problem, getFromMultiApp()->appProblemBase(i));
  }

  _console << "Finished MultiAppAddTransfer " << name() << std::endl;
}
