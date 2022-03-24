//* This file is part of the MOOSE framework
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
