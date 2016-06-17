//Base Classes
#include "FerretApp.h"
#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

//Specific Modules
#include "TensorMechanicsApp.h"
#include "PhaseFieldApp.h"

//Actions
#include "TensorMechanicsActionScaled.h"

//AuxKernels
#include "PolarizationVortexAux.h"
#include "TensorPressureAux.h"
#include "BandGapAuxZnO.h"
#include "BandGapAuxTiO2.h" //should rework this to be a "general" gap kernel
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
#include "ChernSimonsDensityMag.h"
#include "WindingNumberDensity.h"
#include "BandGapAuxZnOwRot.h"

//Boundary Conditions
#include "HydrostaticBC.h"
#include "OpenCircuitBC.h"
#include "StressFreeBC.h"

//Initial Conditions
#include "PerturbedIC.h"
#include "SinIC.h"
#include "FluctuationsIC.h"
#include "AdhocConstIC.h"

//Kernels
#include "ModifiedStressDivergenceTensors.h"
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
#include "FluctuationKernel.h"
#include "FerroelectricCouplingP.h"
//#include "FerroelectricCouplingQ.h"
#include "FerroelectricCouplingX.h"
#include "StressDivergenceTensorsScaled.h"

//Materials
#include "ComputeElectrostrictiveTensor.h"


//Postprocessors
#include "WallEnergy.h"
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

  // ModulesApp::registerObjects(_factory);
  // ModulesApp::associateSyntax(_syntax, _action_factory);

  TensorMechanicsApp::registerObjects(_factory);
  TensorMechanicsApp::associateSyntax(_syntax, _action_factory);

  PhaseFieldApp::registerObjects(_factory);
  PhaseFieldApp::associateSyntax(_syntax, _action_factory);

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
  registerBoundaryCondition(SurfaceMechanicsBC);
  registerBoundaryCondition(HydrostaticBC);
  registerBoundaryCondition(OpenCircuitBC);
  registerBoundaryCondition(StressFreeBC);

  //AuxKernels:
  registerAux(PolarizationVortexAux);
  registerAux(TensorPressureAux);
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
  registerAux(ChernSimonsDensityMag);
  registerAux(CurlP);
  registerAux(BandGapAuxZnO);
  registerAux(BandGapAuxTiO2);
  registerAux(SurfaceChargeAux);
  registerAux(WindingNumberDensity);
  registerAux(BandGapAuxZnOwRot);

  //Kernels

  registerKernel(ModifiedStressDivergenceTensors);
  registerKernel(BulkEnergyDerivativeSixth);
  registerKernel(BulkEnergyDerivativeFourth);
  registerKernel(BulkEnergyDerivativeSixthCoupledT);
  registerKernel(BulkEnergyDerivativeFourthCoupledT);
  registerKernel(WallEnergyDerivative);
  registerKernel(TimeDerivativeScaled);
  registerKernel(FerroelectricCouplingP);
  registerKernel(FluctuationKernel);

//  registerKernel(FerroelectricCouplingQ);
  registerKernel(FerroelectricCouplingX);
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
