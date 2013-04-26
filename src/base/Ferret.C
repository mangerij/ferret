#include "Ferret.h"
#include "FerretApp.h"

#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"
#include "ActionFactory.h"
#include "Parser.h"

#include "Elk.h"

#include "PolarizationVortexAux.h"
#include "PolarizationVortexAuxAction.h"
#include "StressBC.h"
#include "StressFunctionBC.h"
#include "HydrostaticBC.h"
#include "HydrostaticDirichletBC.h"
#include "PolarizationSurfaceCharge.h"
//#include "VortexSurfaceEnergy.h"

namespace Ferret
{
  void registerApps()
  {
    registerApp(FerretApp);
  }

  void registerObjects(Factory & factory)
  {
    registerBoundaryCondition(StressBC);
    registerBoundaryCondition(StressFunctionBC);
    registerBoundaryCondition(HydrostaticBC);
    registerBoundaryCondition(HydrostaticDirichletBC);
    registerBoundaryCondition(PolarizationSurfaceCharge);

    registerAux(PolarizationVortexAux);
    //registerPostprocessor(VortexSurfaceEnergy);
  }
  
  void associateSyntax(Syntax& syntax, ActionFactory & action_factory)
  {
    syntax.registerActionSyntax("PolarizationVortexAuxAction","PolarizationVortexAux");
    registerAction(PolarizationVortexAuxAction, "add_kernel");
  }
}
