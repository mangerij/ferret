//Base Classes
#include "FerretApp.h"
#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

//Specific Modules
#include "TensorMechanicsApp.h"
#include "PhaseFieldApp.h"
#include "MiscApp.h"

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
#include "CurlPMag.h"
#include "BulkEnergyDensity.h"
#include "WallEnergyDensity.h"
#include "SurfaceChargeAux.h"
#include "ChernSimonsDensity.h"
#include "ChernSimonsDensityMag.h"
#include "WindingNumberDensity.h"
#include "BandGapAuxZnOwRot.h"
#include "AngleAux.h"
#include "MieElecFieldReals.h"
#include "MieElecFieldImag.h"
#include "Intensity.h"
#include "OldVar.h"
#include "DielectricTensor.h"


//Boundary Conditions
#include "HydrostaticBC.h"
#include "ScreenedBC.h"
#include "StressFreeBC.h"
#include "StressBC.h"
//#include "ExtrapBC.h"
//#include "MatchedGradValueBC.h"

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
#include "RotatedWallEnergyDerivative.h"
#include "BulkEnergyDerivativeSixth.h"
#include "RotatedBulkEnergyDerivativeSixth.h"
#include "RotatedBulkEnergyDerivative.h"
#include "BulkEnergyDerivativeSixthCoupledT.h"
#include "NoStdBulkEnergyDerivativeSixth.h"
#include "TimeDerivativeScaled.h"
#include "PolarElectricPStrong.h"
#include "PolarElectricEStrong.h"
#include "FluctuationKernel.h"
#include "FerroelectricCouplingP.h"
//#include "FerroelectricCouplingQ.h"
#include "FerroelectricCouplingX.h"
#include "StressDivergenceTensorsScaled.h"
#include "KarmanenkoDriver.h"
#include "KappaTDiffusion.h"
#include "ConstantLatticeMismatch.h"
#include "CoeffParamDiffusion.h"
#include "AnisotropyEnergy.h"
#include "DepolEnergy.h"
#include "SemiconductorChargeCarriers.h"
#include "ThomasFermiPotential.h"
#include "ThomasFermiTerm.h"

//InterfaceKernels
#include "InterfaceDiffusion.h"

//Markers
#include "PolarizationNWEMarker.h"

//Materials
#include "ComputeElectrostrictiveTensor.h"
#include "ComputePhotostrictiveTensor.h"

//Postprocessors
#include "WallEnergy.h"
#include "ThermalEnergy.h"
#include "TotalEnergy.h"
#include "TotalEnergyFlow.h"
#include "TotalEnergyFlowNoElast.h"
#include "TotalEnergyFlowNoElastNoElec.h"
#include "BulkEnergy.h"
#include "BulkEnergyCoupledT.h"
#include "ElectrostaticEnergy.h"
#include "ElasticEnergy.h"
#include "CoupledEnergy.h"
#include "CoupledEnergyCheckShear.h"
#include "GrainSize.h"
#include "DepolarizationEnergy.h"

template<>
InputParameters validParams<FerretApp>()
{
  InputParameters params = validParams<MooseApp>();
  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;

  return params;
}

FerretApp::FerretApp(const InputParameters & parameters) :
    MooseApp(parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  Moose::associateSyntax(_syntax, _action_factory);

  // ModulesApp::registerObjects(_factory);  //uncomment this to activate all modules
  // ModulesApp::associateSyntax(_syntax, _action_factory);

  TensorMechanicsApp::registerObjects(_factory);
  TensorMechanicsApp::associateSyntax(_syntax, _action_factory);

  PhaseFieldApp::registerObjects(_factory);
  PhaseFieldApp::associateSyntax(_syntax, _action_factory);

  MiscApp::registerObjects(_factory);
  MiscApp::associateSyntax(_syntax, _action_factory);

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

  ///BoundaryConditions
  registerBoundaryCondition(SurfaceMechanicsBC);
  registerBoundaryCondition(HydrostaticBC);
  registerBoundaryCondition(ScreenedBC);
  registerBoundaryCondition(StressFreeBC);
  registerBoundaryCondition(StressBC);
  //registerBoundaryCondition(ExtrapBC);
 // registerBoundaryCondition(MatchedGradValueBC);

  ///AuxKernels:
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
  registerAux(CurlPMag);
  registerAux(BandGapAuxZnO);
  registerAux(BandGapAuxTiO2);
  registerAux(SurfaceChargeAux);
  registerAux(WindingNumberDensity);
  registerAux(BandGapAuxZnOwRot);
  registerAux(MieElecFieldReals);
  registerAux(MieElecFieldImag);
  registerAux(Intensity);
  registerAux(AngleAux);
  registerAux(OldVar);
  registerAux(DielectricTensor);

  ///Kernels
  registerKernel(ModifiedStressDivergenceTensors);
  registerKernel(BulkEnergyDerivativeSixth);
  registerKernel(RotatedBulkEnergyDerivativeSixth);
  registerKernel(RotatedBulkEnergyDerivative);
  registerKernel(NoStdBulkEnergyDerivativeSixth);
  registerKernel(BulkEnergyDerivativeSixthCoupledT);
  registerKernel(WallEnergyDerivative);
  registerKernel(RotatedWallEnergyDerivative);
  registerKernel(TimeDerivativeScaled);
  registerKernel(FerroelectricCouplingP);
  registerKernel(FluctuationKernel);
  registerKernel(KarmanenkoDriver);
  registerKernel(KappaTDiffusion);
  registerKernel(ConstantLatticeMismatch);
  registerKernel(AnisotropyEnergy);
  registerKernel(DepolEnergy);
  registerKernel(ThomasFermiPotential);
  registerKernel(ThomasFermiTerm);


  /// registerKernel(FerroelectricCouplingQ);
  registerKernel(FerroelectricCouplingX);
  registerKernel(StressDivergenceTensorsScaled);
  registerKernel(PolarElectricEStrong);
  registerKernel(PolarElectricPStrong);
  registerKernel(Electrostatics);
  registerKernel(CoeffParamDiffusion);
  registerKernel(SemiconductorChargeCarriers);

  ///registerInterfaceKernels
  registerInterfaceKernel(InterfaceDiffusion);

  ///Postprocessors
  registerPostprocessor(BulkEnergy);
  registerPostprocessor(BulkEnergyCoupledT);
  registerPostprocessor(WallEnergy);
  ///registerPostprocessor(ChernSimonsNumber);
  registerPostprocessor(ElectrostaticEnergy);
  registerPostprocessor(TotalEnergy);
  registerPostprocessor(TotalEnergyFlow);
  registerPostprocessor(TotalEnergyFlowNoElast);
  registerPostprocessor(TotalEnergyFlowNoElastNoElec);
  registerPostprocessor(ElasticEnergy);
  registerPostprocessor(CoupledEnergy);
  registerPostprocessor(ThermalEnergy);
  registerPostprocessor(CoupledEnergyCheckShear);
  registerPostprocessor(GrainSize);
  registerPostprocessor(DepolarizationEnergy);

  //Markers
  registerMarker(PolarizationNWEMarker);

  ///Materials
  ///registerMaterial(LinearFerroelectricMaterial); //deprecated
  registerMaterial(ComputeElectrostrictiveTensor);
  registerMaterial(ComputePhotostrictiveTensor);

  ///InitialConditions
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
  registerAction(TensorMechanicsActionScaled, "add_kernel"); ///this is deprecated in our code

#undef registerAction
#define registerAction(tplt, action) action_factory.regLegacy<tplt>(stringifyName(tplt), action)

}
