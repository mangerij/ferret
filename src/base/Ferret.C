#include "Ferret.h"
#include "FerretApp.h"

#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"
#include "ActionFactory.h"
#include "Parser.h"

#include "Elk.h"

#include "PolarizationVortex.h"
#include "PolarizationVortexAction.h"
#include "StressBC.h"
#include "StressFunctionBC.h"
#include "HydrostaticBC.h"
#include "HydrostaticDirichletBC.h"
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
    registerKernel(PolarizationVortex);
    //registerPostprocessor(VortexSurfaceEnergy);
  }
  
  void associateSyntax(Syntax& syntax, ActionFactory & action_factory)
  {
    syntax.registerActionSyntax("PolarizationVortexAction","PolarizationVortex");
    registerAction(PolarizationVortexAction, "add_kernel");
  }
}
