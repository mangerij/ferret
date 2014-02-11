#include "FerretApp.h"
#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"

#include "ElkApp.h"

#include "PolarizationVortexAux.h"
#include "TensorPressureAux.h"
#include "PolarizationVortexAuxAction.h"
#include "StressBC.h"
#include "StressFunctionBC.h"
#include "HydrostaticBC.h"
#include "HydrostaticDirichletBC.h"
#include "SurfaceMechanicsBC.h"
#include "PolarizationSurfaceCharge.h"
//#include "VortexSurfaceEnergy.h"

#include "BulkEnergyDerivative.h"
#include "BulkEnergyDensity.h"
#include "WallEnergyDensity.h"
#include "WallEnergyDerivative.h"
#include "PolarElectricP.h"
#include "PolarElectricPStrong.h"
#include "PolarElectricE.h"
#include "PolarElectricEStrong.h"
#include "BulkEnergy.h"
#include "WallEnergy.h"
#include "ElectricEnergy.h"
#include "ElectrostaticEnergy.h"
#include "TotalEnergy.h"
#include "PolarMaterial.h"
#include "ElectricStatics.h"
#include "PerturbedIC.h"
#include "SinIC.h"
#include "AdhocConstIC.h"
#include "SinFunc.h"
#include "RandomFunc.h"
#include "SphereIC.h"
#include "SphereToCartFunc.h"
#include "ElectrostaticEnergyDensityE.h"
#include "ElectrostaticEnergyDensity.h"
#include "ElectrostaticEnergyDensityCross.h"
#include "ElectrostaticEnergyDensityTotal.h"
#include "PostprocessorAdaptiveDT.h"
#include "CustomDT.h"

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
  registerBoundaryCondition(SurfaceMechanicsBC);
  registerBoundaryCondition(HydrostaticBC);
  registerBoundaryCondition(HydrostaticDirichletBC);
  registerBoundaryCondition(PolarizationSurfaceCharge);

  registerAux(PolarizationVortexAux);
  registerAux(TensorPressureAux);
  registerAux(ElectrostaticEnergyDensity);
  registerAux(ElectrostaticEnergyDensityE);
  registerAux(ElectrostaticEnergyDensityCross);
  registerAux(ElectrostaticEnergyDensityTotal);
  registerAux(BulkEnergyDensity);
  registerAux(WallEnergyDensity);
  //registerPostprocessor(VortexSurfaceEnergy);
  registerKernel(BulkEnergyDerivative);
  registerKernel(WallEnergyDerivative);
  registerKernel(PolarElectricE);
  registerKernel(PolarElectricEStrong);
  registerKernel(PolarElectricP);
  registerKernel(PolarElectricPStrong);
  registerKernel(ElectricStatics);
  registerPostprocessor(BulkEnergy);
  registerPostprocessor(WallEnergy);
  registerPostprocessor(ElectricEnergy);
  registerPostprocessor(ElectrostaticEnergy);
  registerPostprocessor(TotalEnergy);
  registerMaterial(PolarMaterial);
  registerInitialCondition(PerturbedIC);
  registerInitialCondition(SinIC);
  registerInitialCondition(AdhocConstIC);
  registerFunction(SinFunc);
  registerFunction(RandomFunc);
  registerFunction(SphereIC);
  registerFunction(SphereToCartFunc);
  registerTimeStepper(PostprocessorAdaptiveDT);
  registerTimeStepper(CustomDT);
}

void
FerretApp::associateSyntax(Syntax& syntax, ActionFactory & action_factory)
{
  syntax.registerActionSyntax("PolarizationVortexAuxAction","PolarizationVortexAux");
  registerAction(PolarizationVortexAuxAction, "add_kernel");
}
