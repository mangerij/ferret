/**
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <mangeri@fzu.cz>
   and be sure to track new changes at bitbucket.org/mesoscience/ferret

**/

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

//AuxKernels
#include "TensorPressureAux.h"
#include "BandGapAuxZnO.h"
#include "BandGapAuxTiO2.h" //should rework these to be a "general" gap kernel
#include "ElecFieldAux.h"
#include "CurlP.h"
#include "ZCompCurlP.h"
#include "CurlPMag.h"
#include "BulkEnergyDensity.h"
#include "ElectrostrictiveEnergyDensity.h"
#include "ElectrostrictiveCouplingEnergyDensity.h"
#include "RotostrictiveCouplingEnergyDensity.h"
#include "RotoPolarCouplingEnergyDensity.h"
#include "RotoBulkEnergyDensity.h"
#include "WallEnergyDensity.h"
#include "ChernSimonsDensity.h"
#include "ChernSimonsDensityMag.h"
#include "WindingNumberDensity.h"
#include "BandGapAuxZnOwRot.h"
#include "MieElecFieldReals.h"
#include "MieElecFieldImag.h"
#include "Intensity.h"
#include "OldVar.h"
#include "DielectricTensor.h"
#include "NormalizedWallEnergyDensity.h"
#include "PolarMag.h"
#include "DivP.h"
#include "PiezoelectricApprox.h"
#include "Birefringence.h"
#include "AFDWallEnergyDensity.h"
#include "RefractiveIndex.h"
#include "RenormalizedStrain.h"
#include "ChangeInRefractiveIndex.h"
#include "ChangeInRefractiveIndexWithPolar.h"
#include "ChangeInRefractiveIndexWithGCoeffPolar.h"
#include "PkNorm.h"
#include "ChangeInRefractiveIndexElectro.h"
#include "ConvertField.h"
#include "PzSq.h"
#include "MagFieldAux.h"
#include "PlaneAux.h"

//Boundary Conditions
#include "HydrostaticBC.h"
#include "StressBC.h"

//Initial Conditions
#include "FluctuationsIC.h"
#include "RandomConstrainedVectorFieldIC.h"

//Functions
#include "DomainFunc.h"

//Kernels

#include "AnisotropyEnergy.h"

#include "BulkEnergyDerivativeSixth.h"
#include "BulkEnergyDerivativeSixthCoupledT.h"
#include "BulkEnergyDerivativeEighth.h"

#include "BulkAntiferrodistortEnergyDerivativeSixth.h"

#include "TimeDerivativeScaled.h"
#include "PolarElectricPStrong.h"
#include "PolarElectricEStrong.h"
#include "FluctuationKernel.h"
#include "FerroelectricCouplingP.h"
#include "FerroelectricCouplingX.h"
#include "KarmanenkoDriver.h"

#include "LBOBulkEnergyDeriv.h"
#include "DepolEnergy.h"
#include "PontryaginDensity.h"
#include "InPlaneP.h"
#include "ConversePiezoelectricStrain.h"
#include "PiezoelectricStrainCharge.h"
#include "SpinFlexoelectricPolarDerivative.h"
#include "SpinFlexoelectricMagDerivative.h"
#include "SpinFlexoelectricMagGradDerivative.h"
#include "DzyaloshinskiiMagDerivative.h"
#include "DzyaloshinskiiDistortDerivative.h"
#include "MagHStrong.h"
#include "MagMStrong.h"
#include "DampedExchangeLLG.h"
#include "DampedAnisotropyLLG.h"
#include "DampedConstrainedExchangeLLG.h"
#include "DampedConstrainedAnisotropyLLG.h"
#include "LagrangeMagConstraint.h"
#include "LagrangeLambdaConstraint.h"
#include "LocalSaturationConstraint.h"
#include "AFDAntiphaseEnergyDerivative.h"

#include "SurfaceMechanicsBC.h" //not sure why this is called a BC
#include "Electrostatics.h"

#include "WallEnergyDerivative.h"
#include "WallEnergyDerivativeAlt.h"

#include "RotostrictiveCouplingDistortDerivative.h"
#include "RotostrictiveCouplingDispDerivative.h"
#include "ElectrostrictiveCouplingPolarDerivative.h"
#include "ElectrostrictiveCouplingDispDerivative.h"
#include "RotoPolarCoupledEnergyPolarDerivativeAlt.h"
#include "RotoBulkEnergyDerivativeEighthAlt.h"
#include "RotoPolarCoupledEnergyDistortDerivativeAlt.h"
#include "ConstField.h"

//InterfaceKernels
#include "InterfaceDiffusion.h"

//Markers
#include "PolarizationNWEMarker.h"

//Materials
#include "ComputeElectrostrictiveTensor.h"
#include "ComputeElastoopticTensor.h"
#include "ComputeDeltaIndicatrix.h"
#include "ComputeIndicatrix.h"
#include "ComputePolarOpticTensor.h"
#include "ComputePolarOpticGCoeffTensor.h"
#include "ComputeElectroopticTensor.h"
#include "ComputeGCoeffTensor.h"
#include "ComputeDeltaIndicatrixElectro.h"
#include "ComputePiezoTensor.h"
#include "ComputePiezostrictiveTensor.h"
#include "ComputeRankTwoLandauTensor.h"
#include "ComputeRankFourLandauTensor.h"

//Postprocessors
#include "WallEnergy.h"
#include "ThermalEnergy.h"

#include "PostprocessorsSum.h"

#include "BulkEnergy.h"
#include "BulkEnergyCoupledT.h"
#include "ElectrostaticEnergy.h"
#include "ExtElectrostaticEnergy.h"
#include "ElasticEnergy.h"
#include "ElectrostrictiveEnergy.h"
#include "DepolarizationEnergy.h"
#include "AnisotropicEnergy.h"

#include "TotalWinding.h"
#include "EnergyRatePostprocessor.h"
#include "LBOBulkEnergy.h"
#include "BulkAntiferrodistortEnergy.h"
#include "MagneticExchangeEnergy.h"
#include "MagneticAnisotropyEnergy.h"
#include "DMInteractionEnergy.h"
#include "RotopolarCouplingEnergy.h"
#include "BulkEnergyEighth.h"
#include "RotoBulkEnergyEighth.h"
#include "RotoPolarCoupledEnergyEighth.h"
#include "AFDWallEnergy.h"
#include "RotostrictiveCouplingEnergy.h"
#include "ElectrostrictiveCouplingEnergy.h"
#include "PolarizationValue.h"
#include "PolarizationComponentValue.h"

//UserObjects

//#include "GrainBoundaryUserObject.h"

//VectorPostprocessors

//#include "GrainBoundaryOutwardNormalLineValueSampler.h"


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
  registerBoundaryCondition(StressBC);

  ///AuxKernels:
  registerAux(TensorPressureAux);
  registerAux(BulkEnergyDensity);
  registerAux(ElectrostrictiveEnergyDensity);
  registerAux(ElectrostrictiveCouplingEnergyDensity);
  registerAux(RotostrictiveCouplingEnergyDensity);
  registerAux(RotoPolarCouplingEnergyDensity);
  registerAux(RotoBulkEnergyDensity);
  registerAux(WallEnergyDensity);
  registerAux(AFDWallEnergyDensity);
  registerAux(ElecFieldAux);
  registerAux(ChernSimonsDensity);
  registerAux(ChernSimonsDensityMag);
  registerAux(CurlP);
  registerAux(ZCompCurlP);
  registerAux(CurlPMag);
  registerAux(BandGapAuxZnO);
  registerAux(BandGapAuxTiO2);
  registerAux(WindingNumberDensity);
  registerAux(BandGapAuxZnOwRot);
  registerAux(MieElecFieldReals);
  registerAux(MieElecFieldImag);
  registerAux(Intensity);
  registerAux(OldVar);
  registerAux(DielectricTensor);
  registerAux(NormalizedWallEnergyDensity);
  registerAux(PolarMag);
  registerAux(DivP);
  registerAux(PiezoelectricApprox);
  registerAux(Birefringence);
  registerAux(ChangeInRefractiveIndex);
  registerAux(RefractiveIndex);
  registerAux(RenormalizedStrain);
  registerAux(ChangeInRefractiveIndexWithPolar);
  registerAux(ChangeInRefractiveIndexWithGCoeffPolar);
  registerAux(PkNorm);
  registerAux(PontryaginDensity);
  registerAux(MagFieldAux);
  registerAux(InPlaneP);
  registerAux(ChangeInRefractiveIndexElectro);
  registerAux(ConvertField);
  registerAux(PzSq);
  registerAux(PlaneAux);

  registerFunction(DomainFunc);

  ///Kernels
  registerKernel(BulkEnergyDerivativeSixth);
  registerKernel(BulkEnergyDerivativeSixthCoupledT);
  registerKernel(WallEnergyDerivative);
  registerKernel(WallEnergyDerivativeAlt);
  registerKernel(TimeDerivativeScaled);
  registerKernel(FerroelectricCouplingP);
  registerKernel(FluctuationKernel);
  registerKernel(KarmanenkoDriver);
  registerKernel(AnisotropyEnergy);
  registerKernel(LBOBulkEnergyDeriv);
  registerKernel(DepolEnergy);
  registerKernel(FerroelectricCouplingX);
  registerKernel(PolarElectricEStrong);
  registerKernel(PolarElectricPStrong);
  registerKernel(Electrostatics);
  registerKernel(ConversePiezoelectricStrain);
  registerKernel(PiezoelectricStrainCharge);
  //experimental magnetic and magnetoelectric kernels under development
  // see Popkov et al PRB,92,140414(R) (2015):
  registerKernel(SpinFlexoelectricPolarDerivative);
  registerKernel(SpinFlexoelectricMagDerivative);
  registerKernel(SpinFlexoelectricMagGradDerivative);
  registerKernel(DzyaloshinskiiMagDerivative);
  registerKernel(DzyaloshinskiiDistortDerivative);
  registerKernel(MagMStrong);
  registerKernel(MagHStrong);
  registerKernel(BulkAntiferrodistortEnergyDerivativeSixth);
  registerKernel(RotoPolarCoupledEnergyPolarDerivativeAlt);
  registerKernel(AFDAntiphaseEnergyDerivative);
  registerKernel(BulkEnergyDerivativeEighth);
  registerKernel(RotostrictiveCouplingDistortDerivative);
  registerKernel(RotostrictiveCouplingDispDerivative);
  registerKernel(ElectrostrictiveCouplingPolarDerivative);
  registerKernel(ElectrostrictiveCouplingDispDerivative);
  registerKernel(RotoBulkEnergyDerivativeEighthAlt);
  registerKernel(RotoPolarCoupledEnergyDistortDerivativeAlt);
  registerKernel(ConstField);
  registerKernel(DampedAnisotropyLLG);
  registerKernel(DampedExchangeLLG);
  registerKernel(DampedConstrainedAnisotropyLLG);
  registerKernel(DampedConstrainedExchangeLLG);
  registerKernel(LagrangeMagConstraint);
  registerKernel(LagrangeLambdaConstraint);
  registerKernel(LocalSaturationConstraint);

  ///registerInterfaceKernels
  registerInterfaceKernel(InterfaceDiffusion); //used for testing purposes

  ///Postprocessors
  registerPostprocessor(BulkEnergy);
  registerPostprocessor(BulkEnergyCoupledT);
  registerPostprocessor(WallEnergy);
  registerPostprocessor(ElectrostaticEnergy);
  registerPostprocessor(ElasticEnergy);
  registerPostprocessor(ElectrostrictiveEnergy); //used to be called CoupledEnergy
  registerPostprocessor(ThermalEnergy);
  registerPostprocessor(DepolarizationEnergy);
  registerPostprocessor(AnisotropicEnergy);
  registerPostprocessor(TotalWinding);

  registerPostprocessor(PostprocessorsSum);

  registerPostprocessor(EnergyRatePostprocessor);
  registerPostprocessor(LBOBulkEnergy);
  registerPostprocessor(ExtElectrostaticEnergy);
  registerPostprocessor(BulkAntiferrodistortEnergy);
  registerPostprocessor(MagneticExchangeEnergy);
  registerPostprocessor(MagneticAnisotropyEnergy);
  registerPostprocessor(DMInteractionEnergy);
  registerPostprocessor(RotopolarCouplingEnergy);
  registerPostprocessor(BulkEnergyEighth);
  registerPostprocessor(RotoBulkEnergyEighth);
  registerPostprocessor(RotoPolarCoupledEnergyEighth);
  registerPostprocessor(AFDWallEnergy);
  registerPostprocessor(RotostrictiveCouplingEnergy);
  registerPostprocessor(ElectrostrictiveCouplingEnergy);
  registerPostprocessor(PolarizationValue);
  registerPostprocessor(PolarizationComponentValue);

   //VectorPostprocessors
 // registerVectorPostprocessor(GrainBoundaryOutwardNormalLineValueSampler);
  //UserObjects
  
  //registerUserObject(GrainBoundaryUserObject);

  //Markers
  registerMarker(PolarizationNWEMarker);

  ///Materials
  ///registerMaterial(LinearFerroelectricMaterial); //deprecated, long live this simple method!
  registerMaterial(ComputeElectrostrictiveTensor);
  registerMaterial(ComputeElastoopticTensor);
  registerMaterial(ComputeDeltaIndicatrix);
  registerMaterial(ComputeIndicatrix);
  registerMaterial(ComputePolarOpticTensor);
  registerMaterial(ComputePolarOpticGCoeffTensor);
  registerMaterial(ComputeElectroopticTensor);
  registerMaterial(ComputeGCoeffTensor);
  registerMaterial(ComputeDeltaIndicatrixElectro);
  registerMaterial(ComputePiezoTensor);
  registerMaterial(ComputePiezostrictiveTensor);
  registerMaterial(ComputeRankTwoLandauTensor);
  registerMaterial(ComputeRankFourLandauTensor);


  ///InitialConditions
  registerInitialCondition(FluctuationsIC);
  registerInitialCondition(RandomConstrainedVectorFieldIC);


#undef registerObject
#define registerObject(name) factory.regLegacy<name>(stringifyName(name))
}
