//Base Classes
#include "FerretApp.h"
#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"
#include "ModulesApp.h"

//Actions
#include "TensorMechanicsActionScaled.h"
#include "PolarizationVortexAuxAction.h"

//AuxKernels
#include "PolarizationVortexAux.h"
#include "TensorPressureAux.h"
#include "BandGapAuxZnO.h"
#include "BandGapAuxTiO2.h"
#include "Ex_fieldAux.h"
#include "Ey_fieldAux.h"
#include "Ez_fieldAux.h"
#include "Px_fieldAux.h"
#include "Py_fieldAux.h"
#include "Pz_fieldAux.h"
#include "BoundCharge.h"
#include "BulkEnergyDensity.h"
#include "WallEnergyDensity.h"
#include "SurfaceChargeAux.h"
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

//Initial Conditions
#include "PerturbedIC.h"
#include "SinIC.h"
#include "AdhocConstIC.h"

//Kernels
#include "SurfaceMechanicsBC.h" //not sure why this is called a BC
#include "Electrostatics.h"
#include "WallEnergyDerivative.h"
#include "WallEnergyDerivative_scaled.h" //deprecated
#include "BulkEnergyDerivative.h"
#include "BulkEnergyDerivative_nosixth.h"
#include "BulkEnergyDerivative_scaled.h" //deprecated
#include "TimeDerivativeScaled.h"
#include "PolarElectricP.h"
#include "PolarElectricPStrong.h"
#include "PolarElectricE.h"
#include "PolarElectricEStrong.h"
#include "FerroelectricCouplingP.h"
#include "FerroelectricCouplingU.h"
#include "StressDivergenceTensorsScaled.h"

//Materials
#include "LinearFerroelectricMaterial.h"
#include "PolarMaterial.h"

//Postprocessors
#include "WallEnergy.h"
#include "TotalEnergy.h"
#include "BulkEnergy.h"
#include "ElectricEnergy.h" //deprecated
#include "ElectrostaticEnergy.h"
#include "TotalEnergyGradient.h" //deprecated
#include "TotalEnergyGradientL2.h"
#include "PercentChangePostprocessor.h"//added to MOOSE 6/23, so deprecated here
#include "ElasticEnergy.h"

//Time steppers
#include "PostprocessorAdaptiveDT.h"
#include "CustomDT.h"
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
  registerAux(Ex_fieldAux);
  registerAux(Ey_fieldAux);
  registerAux(Ez_fieldAux);
  registerAux(Px_fieldAux);
  registerAux(Py_fieldAux);
  registerAux(Pz_fieldAux);
  registerAux(BoundCharge);
  registerAux(BandGapAuxZnO);
  registerAux(BandGapAuxTiO2);

  registerAux(SurfaceChargeAux);
  //registerPostprocessor(VortexSurfaceEnergy);
  registerKernel(BulkEnergyDerivative);
  registerKernel(BulkEnergyDerivative_nosixth);
  registerKernel(WallEnergyDerivative);

  registerKernel(BulkEnergyDerivative_scaled);
  registerKernel(WallEnergyDerivative_scaled);
  registerKernel(TimeDerivativeScaled);
  registerKernel(FerroelectricCouplingP);
  registerKernel(FerroelectricCouplingU);
  registerKernel(StressDivergenceTensorsScaled);

  registerKernel(PolarElectricE);
  registerKernel(PolarElectricEStrong);
  registerKernel(PolarElectricP);
  registerKernel(PolarElectricPStrong);
  registerKernel(Electrostatics);

  registerPostprocessor(BulkEnergy);
  registerPostprocessor(WallEnergy);
  registerPostprocessor(ElectricEnergy);
  registerPostprocessor(ElectrostaticEnergy);
  registerPostprocessor(TotalEnergy);
  registerPostprocessor(TotalEnergyGradient);
  registerPostprocessor(TotalEnergyGradientL2);
  registerPostprocessor(PercentChangePostprocessor);
  registerPostprocessor(ElasticEnergy);

  registerMaterial(PolarMaterial);
  registerMaterial(LinearFerroelectricMaterial);

  registerInitialCondition(PerturbedIC);
  registerInitialCondition(SinIC);
  registerInitialCondition(AdhocConstIC);

  registerFunction(SinFunc);
  registerFunction(RandomFunc);
  registerFunction(SphereIC);
  registerFunction(SphereToCartFunc);

  registerTimeStepper(PostprocessorAdaptiveDT);
  registerTimeStepper(CustomDT);
  registerTimeStepper(TransientHalf);
}

void
FerretApp::associateSyntax(Syntax& syntax, ActionFactory & action_factory)
{
  syntax.registerActionSyntax("PolarizationVortexAuxAction","PolarizationVortexAux");
  registerAction(PolarizationVortexAuxAction, "add_kernel");

  syntax.registerActionSyntax("TensorMechanicsActionScaled", "Kernels/TensorMechanicsScaled");
  registerAction(TensorMechanicsActionScaled, "add_kernel");

}
