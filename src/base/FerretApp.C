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

#include "AppFactory.h"
#include "FerretApp.h"
#include "Moose.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

// Optional dependence on AnotherApp (maybe ScalFMM?)
//#ifdef ANOTHER_ENABLED
//#include "AnotherApp.h"
//#endif

InputParameters
FerretApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("error_unused") = false;
  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

// When using the new Registry system, this line is required so that
// dependent apps know about the MastodonApp label.
registerKnownLabel("FerretApp");

FerretApp::FerretApp(InputParameters parameters) : MooseApp(parameters)
{
  FerretApp::registerAll(_factory, _action_factory, _syntax);
}

void
FerretApp::registerApps()
{
  registerApp(FerretApp);
}

void
FerretApp::registerAll(Factory & factory, ActionFactory & action_factory, Syntax & syntax)
{

  Registry::registerObjectsTo(factory, {"FerretApp"});
  Registry::registerActionsTo(action_factory, {"FerretApp"});
  ModulesApp::registerAll(factory, action_factory, syntax);
//#ifdef ANOTHER_ENABLED
//  AnotherApp::registerAll(factory, action_factory, syntax);
//#endif

  syntax.registerActionSyntax("ABO3CoupledPhaseFieldAction", "Ferret/ABO3CoupledPhaseField");

}

extern "C" void
FerretApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  FerretApp::registerAll(f, af, s);
}
