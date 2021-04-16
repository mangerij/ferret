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

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "LKScalarAction.h"
#include "Factory.h"
#include "Parser.h"
#include "FEProblem.h"

#include "AddVariableAction.h"

#include "libmesh/string_to_enum.h"


registerMooseAction("FerretApp", LKScalarAction, "add_scalar_kernel");


template <>
InputParameters
validParams<LKScalarAction>()
{
  InputParameters params = validParams<Action>();
  params.addClassDescription("Set up Landau-Khalatnikov ScalarKernels");
  params.addRequiredParam<std::vector<NonlinearVariableName>>(
      "variables", "The names of order parameter and strain tensor components");
  return params;
}

LKScalarAction::LKScalarAction(InputParameters params) : Action(params)
{
}

void
LKScalarAction::act()
{
  if (_current_task == "add_scalar_kernel")
  {
    unsigned int _ord_num = 12;
    for (unsigned int kk = 0; kk < _ord_num; ++kk)
    {
      if(_polar_time_dependence==true)
      {
        InputParameters params = _factory.getValidParams("ODETimeDerivative");
        params.set<NonlinearVariableName>("variable") = variables[kk];
        params.applyParameters(parameters());

        std::string s_kernel_name = "odetdp_" + Moose::stringify(kk);
        _problem->addKernel("ODETimeDerivative", s_kernel_name, params);
      }
    }
  }
}

