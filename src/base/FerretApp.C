#include "FerretApp.h"
#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"

#include "ElkApp.h"

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
#include "SinIC.h"
#include "AdhocConstIC.h"
#include "SinFunc.h"
#include "RandomFunc.h"
#include "SphereIC.h"
#include "SphereToCartFunc.h"

template<>
InputParameters validParams<FerretApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

FerretApp::FerretApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(libMesh::processor_id());
  
  Moose::registerObjects(_factory);
  Moose::associateSyntax(_syntax, _action_factory);
  
  ElkApp::registerObjects(_factory);
  ElkApp::associateSyntax(_syntax, _action_factory);
  
  FerretApp::registerObjects(_factory);
  FerretApp::associateSyntax(_syntax, _action_factory);
}

FerretApp::~FerretApp()
{
}

void
FerretApp::registerApps()
{
  registerApp(FerretApp);
}

void
FerretApp::registerObjects(Factory & factory)
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
  registerInitialCondition(SinIC);
  registerInitialCondition(AdhocConstIC);
  registerFunction(SinFunc);
  registerFunction(RandomFunc);
  registerFunction(SphereIC);
  registerFunction(SphereToCartFunc);
}

void
FerretApp::associateSyntax(Syntax& syntax, ActionFactory & action_factory)
{
  syntax.registerActionSyntax("PolarizationVortexAuxAction","PolarizationVortexAux");
  registerAction(PolarizationVortexAuxAction, "add_kernel");
}
