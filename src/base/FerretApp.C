//Base Classes
#include "FerretApp.h"
#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

//Actions
#include "TensorMechanicsActionScaled.h"

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
#include "CurlP.h"
#include "BulkEnergyDensity.h"
#include "WallEnergyDensity.h"
#include "SurfaceChargeAux.h"
#include "ChernSimonsDensity.h"
#include "ScreenAux.h"
#include "PolarizationSurfaceCharge.h" // What was this?
#include "WindingNumberDensity.h"

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
#include "OpenCircuitBC.h"


//Initial Conditions
#include "PerturbedIC.h"
#include "SinIC.h"
#include "FluctuationsIC.h"
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
//#include "FerroelectricCouplingU.h" removed, but we can leave this here as a testament to mediocrity
//#include "FerroelectricCouplingX.h"
#include "StressDivergenceTensorsScaled.h"

//Materials
//#include "LinearFerroelectricMaterial.h"
//new "plug and play" approach for electrostriction
#include "ComputeElectrostrictiveTensor.h"
//#include "ComputeRotatedElectrostrictiveTensorBase.h"
//#include "ComputeElectrostrictiveTensorBase.h"
//#include "PolarMaterial.h"

//Postprocessors
#include "WallEnergy.h"
#include "PZTWallEnergy.h"
//#include "WallEnergyFourth.h"
#include "TotalEnergy.h"
#include "TotalEnergyFlow.h"
#include "BulkEnergy.h"
#include "BulkEnergyFourth.h"
#include "ElectrostaticEnergy.h"
#include "ElasticEnergy.h"
#include "CoupledEnergy.h"

template<>
InputParameters validParams<FerretApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

FerretApp::FerretApp(const InputParameters & parameters) :
    MooseApp(parameters)
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

// void
//FerretApp::registerApps()
//{
//  registerApp(FerretApp);
//}

void
FerretApp::registerApps()
{
#undef  registerApp
#define registerApp(name) AppFactory::instance().reg<name>(#name)
  registerApp(FerretApp);
#undef  registerApp
#define registerApp(name) AppFactory::instance().regLegacy<name>(#name)
}



void
FerretApp::registerObjects(Factory & factory)
{

#undef registerObject
#define registerObject(name) factory.reg<name>(stringifyName(name))

  //BoundaryConditions
  registerBoundaryCondition(StressBC);
  registerBoundaryCondition(StressFunctionBC);
  registerBoundaryCondition(SurfaceMechanicsBC);
  registerBoundaryCondition(HydrostaticBC);
  registerBoundaryCondition(HydrostaticDirichletBC);
  registerBoundaryCondition(PolarizationSurfaceCharge);
  registerBoundaryCondition(DepolScreenBC);
  registerBoundaryCondition(OpenCircuitBC);

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
  registerAux(ChernSimonsDensity);
  registerAux(CurlP);
  registerAux(BandGapAuxZnO);
  registerAux(BandGapAuxTiO2);
  registerAux(SurfaceChargeAux);
  registerAux(ScreenAux);
  registerAux(WindingNumberDensity);

  //Kernels
  registerKernel(BulkEnergyDerivativeSixth);
  registerKernel(BulkEnergyDerivativeFourth);
  registerKernel(BulkEnergyDerivativeSixthCoupledT);
  registerKernel(BulkEnergyDerivativeFourthCoupledT);
  registerKernel(WallEnergyDerivative);
  registerKernel(TimeDerivativeScaled);
  registerKernel(FerroelectricCouplingP);
 // registerKernel(FerroelectricCouplingU);
 // registerKernel(FerroelectricCouplingX);
  registerKernel(StressDivergenceTensorsScaled);
  registerKernel(PolarElectricEStrong);
  registerKernel(PolarElectricPStrong);
  registerKernel(Electrostatics);

  //Postprocessors
  registerPostprocessor(BulkEnergy);
  registerPostprocessor(WallEnergy);
  //registerPostprocessor(ChernSimonsNumber);
  registerPostprocessor(BulkEnergyFourth);
  registerPostprocessor(ElectrostaticEnergy);
  registerPostprocessor(TotalEnergy);
  registerPostprocessor(TotalEnergyFlow);
  registerPostprocessor(ElasticEnergy);
  registerPostprocessor(CoupledEnergy);

  //Materials
  //registerMaterial(LinearFerroelectricMaterial);

  registerMaterial(ComputeElectrostrictiveTensor);
  //registerMaterial(ComputeRotatedElectrostrictiveTensorBase);
  //registerMaterial(ComputeElectrostrictiveTensorBase);

  //InitialConditions
  registerInitialCondition(PerturbedIC);
  registerInitialCondition(SinIC);
  registerInitialCondition(FluctuationsIC);
  registerInitialCondition(AdhocConstIC);

#undef registerObject
#define registerObject(name) factory.regLegacy<name>(stringifyName(name))
}

void
FerretApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{

#undef registerAction
#define registerAction(tplt, action) action_factory.reg<tplt>(stringifyName(tplt), action)

  syntax.registerActionSyntax("TensorMechanicsActionScaled", "Kernels/TensorMechanicsScaled");
  registerAction(TensorMechanicsActionScaled, "add_kernel"); //this is deprecated in our code

#undef registerAction
#define registerAction(tplt, action) action_factory.regLegacy<tplt>(stringifyName(tplt), action)

}
