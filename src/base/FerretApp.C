//Base Classes
#include "FerretApp.h"
#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"
#include "ModulesApp.h"

//Actions
#include "TensorMechanicsActionScaled.h" //not sure the tensormechanics action works right now
#include "PolarizationVortexAuxAction.h"

//AuxKernels
#include "PolarizationVortexAux.h"
#include "TensorPressureAux.h"
#include "BandGapAuxZnO.h"
#include "BandGapAuxTiO2.h"
#include "ExFieldAux.h"
#include "EyFieldAux.h"
#include "EzFieldAux.h"
#include "PxFieldAux.h"
#include "PyFieldAux.h"
#include "PzFieldAux.h"
#include "BoundCharge.h"
#include "BulkEnergyDensity.h"
#include "WallEnergyDensity.h"
#include "SurfaceChargeAux.h"
#include "ScreenAux.h"
#include "PolarizationSurfaceCharge.h" //?

//not sure what these are
#include "ElectrostaticEnergyDensityE.h"
#include "ElectrostaticEnergyDensity.h"
#include "ElectrostaticEnergyDensityCross.h"
#include "ElectrostaticEnergyDensityTotal.h"

//#include "VortexSurfaceEnergy.h"

//Boundary Conditions
#include "StressBC.h"
#include "StressFunctionBC.h"
#include "HydrostaticBC.h"
#include "HydrostaticDirichletBC.h"
#include "DepolScreenBC.h"

//Initial Conditions
#include "PerturbedIC.h"
#include "SinIC.h"
#include "AdhocConstIC.h"

//Kernels
#include "SurfaceMechanicsBC.h" //not sure why this is called a BC
#include "Electrostatics.h"
#include "WallEnergyDerivative.h"
#include "BulkEnergyDerivativeSixth.h"
#include "BulkEnergyDerivativeFourth.h"
#include "BulkEnergyDerivativeSixthCoupledT.h"
#include "BulkEnergyDerivativeFourthCoupledT.h"
#include "TimeDerivativeScaled.h"
#include "PolarElectricPStrong.h"
#include "PolarElectricEStrong.h"
#include "FerroelectricCouplingP.h"
#include "FerroelectricCouplingU.h"
#include "StressDivergenceTensorsScaled.h"

//Materials
#include "LinearFerroelectricMaterial.h"
#include "PolarMaterial.h"

//Postprocessors
#include "WallEnergy.h"
//#include "WallEnergyFourth.h"
#include "TotalEnergy.h"
#include "BulkEnergy.h"
#include "BulkEnergyFourth.h"
#include "ElectrostaticEnergy.h"
#include "TotalEnergyGradient.h" //deprecated
#include "TotalEnergyGradientL2.h"
#include "ElasticEnergy.h"
#include "CoupledEnergy.h"

//Time steppers
#include "PostprocessorAdaptiveDT.h"
#include "TransientHalf.h" //these are all junk

//custom functions
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
  srand(processor_id());

  Moose::registerObjects(_factory);
  Moose::associateSyntax(_syntax, _action_factory);

  ModulesApp::registerObjects(_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);

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
  //BoundaryConditions
  registerBoundaryCondition(StressBC);
  registerBoundaryCondition(StressFunctionBC);
  registerBoundaryCondition(SurfaceMechanicsBC);
  registerBoundaryCondition(HydrostaticBC);
  registerBoundaryCondition(HydrostaticDirichletBC);
  registerBoundaryCondition(PolarizationSurfaceCharge);
  registerBoundaryCondition(DepolScreenBC);

  //AuxKernels:
  registerAux(PolarizationVortexAux);
  registerAux(TensorPressureAux);
  registerAux(ElectrostaticEnergyDensity);
  registerAux(ElectrostaticEnergyDensityE);
  registerAux(ElectrostaticEnergyDensityCross);
  registerAux(ElectrostaticEnergyDensityTotal);
  registerAux(BulkEnergyDensity);
  registerAux(WallEnergyDensity);
  registerAux(ExFieldAux);
  registerAux(EyFieldAux);
  registerAux(EzFieldAux);
  registerAux(PxFieldAux);
  registerAux(PyFieldAux);
  registerAux(PzFieldAux);
  registerAux(BoundCharge);
  registerAux(BandGapAuxZnO);
  registerAux(BandGapAuxTiO2);
  registerAux(SurfaceChargeAux);
  registerAux(ScreenAux);
  //registerPostprocessor(VortexSurfaceEnergy);

  //Kernels
  registerKernel(BulkEnergyDerivativeSixth);
  registerKernel(BulkEnergyDerivativeFourth);
  registerKernel(BulkEnergyDerivativeSixthCoupledT);
  registerKernel(BulkEnergyDerivativeFourthCoupledT);
  registerKernel(WallEnergyDerivative);
  registerKernel(TimeDerivativeScaled);
  registerKernel(FerroelectricCouplingP);
  registerKernel(FerroelectricCouplingU);
  registerKernel(StressDivergenceTensorsScaled);
  registerKernel(PolarElectricEStrong);
  registerKernel(PolarElectricPStrong);
  registerKernel(Electrostatics);

  //Postprocessors
  registerPostprocessor(BulkEnergy);
  registerPostprocessor(WallEnergy);
  registerPostprocessor(BulkEnergyFourth);
  registerPostprocessor(ElectrostaticEnergy);
  registerPostprocessor(TotalEnergy);
  registerPostprocessor(TotalEnergyGradient);
  registerPostprocessor(TotalEnergyGradientL2);
 // registerPostprocessor(PercentChangePostprocessor); //added to MOOSE (deprecated)
  registerPostprocessor(ElasticEnergy);
  registerPostprocessor(CoupledEnergy);

  //Materials
  registerMaterial(PolarMaterial);
  registerMaterial(LinearFerroelectricMaterial);

  //InitialConditions
  registerInitialCondition(PerturbedIC);
  registerInitialCondition(SinIC);
  registerInitialCondition(AdhocConstIC);

  //Functions
  registerFunction(SinFunc);
  registerFunction(RandomFunc);
  registerFunction(SphereIC);
  registerFunction(SphereToCartFunc);

  //TimeStepper
  registerTimeStepper(PostprocessorAdaptiveDT);
  registerTimeStepper(TransientHalf);
}

void
FerretApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  syntax.registerActionSyntax("PolarizationVortexAuxAction","PolarizationVortexAux");
  registerAction(PolarizationVortexAuxAction, "add_kernel");

  syntax.registerActionSyntax("TensorMechanicsActionScaled", "Kernels/TensorMechanicsScaled");
  registerAction(TensorMechanicsActionScaled, "add_kernel");

}
