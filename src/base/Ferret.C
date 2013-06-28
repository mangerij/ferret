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

#include "BulkEnergyDerivative.h"
#include "WallEnergyDerivative.h"
#include "PolarElectricP.h"
#include "PolarElectricE.h"
#include "BulkEnergy.h"
#include "PolarMaterial.h"
#include "ElectricStatics.h"
#include "PerturbedIC.h"
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
    registerKernel(BulkEnergyDerivative);
    registerKernel(WallEnergyDerivative);
    registerKernel(PolarElectricE);
    registerKernel(PolarElectricP);
    registerKernel(ElectricStatics);
    registerPostprocessor(BulkEnergy);
    registerMaterial(PolarMaterial);
    registerInitialCondition(PerturbedIC);
  }

  void associateSyntax(Syntax& syntax, ActionFactory & action_factory)
  {
    syntax.registerActionSyntax("PolarizationVortexAuxAction","PolarizationVortexAux");
    registerAction(PolarizationVortexAuxAction, "add_kernel");
  }
}
