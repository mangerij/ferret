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

   For help with FERRET please contact J. Mangeri <john.mangeri@uconn.edu>
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
#include "PolarizationVortexAux.h"
#include "TensorPressureAux.h"
#include "BandGapAuxZnO.h"
#include "BandGapAuxTiO2.h" //should rework these to be a "general" gap kernel
#include "ElecFieldAux.h"
#include "CurlP.h"
#include "ZCompCurlP.h"
#include "CurlPMag.h"
#include "BulkEnergyDensity.h"
#include "WallEnergyDensity.h"
#include "SurfaceChargeAux.h"
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
#include "SemiconductingChargeCarriersAux.h"
#include "Birefringence.h"
#include "RefractiveIndex.h"
#include "ChangeInRefractiveIndex.h"
#include "ChangeInRefractiveIndexWithPolar.h"
#include "ChangeInRefractiveIndexWithGCoeffPolar.h"
#include "PkNorm.h"
#include "SemiconductingChargeCarriersPolyLogAux.h"
#include "ChangeInRefractiveIndexElectro.h"

//Boundary Conditions
#include "HydrostaticBC.h"
#include "StressBC.h"

//Initial Conditions
#include "FluctuationsIC.h"

//Kernels
#include "ModifiedStressDivergenceTensors.h"
#include "SurfaceMechanicsBC.h" //not sure why this is called a BC
#include "Electrostatics.h"
#include "WallEnergyDerivative.h"
#include "WallEnergyDerivativeAlt.h"
#include "RotatedWallEnergyDerivative.h"
#include "BulkEnergyDerivativeSixth.h"
#include "BulkEnergyDerivativeSixthAlt.h"
#include "BulkEnergyDerivativePSTO.h"
#include "RotatedBulkEnergyDerivativeSixth.h"
#include "RotatedBulkEnergyDerivative.h"
#include "BulkEnergyDerivativeSixthCoupledT.h"
#include "NoStdBulkEnergyDerivativeSixth.h"
#include "TimeDerivativeScaled.h"
#include "PolarElectricPStrong.h"
#include "PolarElectricEStrong.h"
#include "FluctuationKernel.h"
#include "FerroelectricCouplingP.h"
#include "FerroelectricCouplingX.h"
#include "StressDivergenceTensorsScaled.h"
#include "KarmanenkoDriver.h"
#include "KappaTDiffusion.h"
#include "CoeffParamDiffusion.h"
#include "AnisotropyEnergy.h"
#include "LBOBulkEnergyDeriv.h"
#include "DepolEnergy.h"
#include "SemiconductorChargeCarriers.h"
#include "ThomasFermiPotential.h"
#include "ThomasFermiTerm.h"
#include "RenormalizedFreeEnergy.h"
#include "AnisotropicElectrostatics.h"
#include "NerstPlanckDrivingTerm.h"
#include "NerstPlanckDiffusive.h"
#include "FreeChargeContribution.h"
#include "HoleChargeContribution.h"
#include "AcceptorIonContribution.h"
#include "SkyrmionChargeDensityZ.h"
#include "SemiconductorChargeCarriersPolyLog.h"
#include "PolarElectricEStrongAlt.h"
#include "ConversePiezoelectricStrain.h"
#include "PiezoelectricStrainCharge.h"

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

//Postprocessors
#include "WallEnergy.h"
#include "ThermalEnergy.h"
#include "TotalEnergy.h"
#include "TotalEnergyPSTO.h"
#include "TotalEnergyPSTOcoupled.h"
#include "CoupledEnergyPSTO.h"
#include "TotalEnergyFlow.h"
#include "TotalEnergyFlowNoElast.h"
#include "TotalEnergyFlowNoElastNoElec.h"
#include "BulkEnergy.h"
#include "BulkEnergyPSTO.h"
#include "BulkEnergyCoupledT.h"
#include "ElectrostaticEnergy.h"
#include "ExtElectrostaticEnergy.h"
#include "ElasticEnergy.h"
#include "CoupledEnergy.h"
#include "ElectrostrictiveEnergy.h"
#include "CoupledEnergyCheckShear.h"
#include "GrainSize.h"
#include "DepolarizationEnergy.h"
#include "AnisotropicEnergy.h"
#include "TotalEnergyG.h"
#include "RenormalizedBulkEnergy.h"
#include "TotalEnergyP.h"
#include "TotalEnergySkFlow.h"
#include "TotalWinding.h"
#include "EnergyRatePostprocessor.h"
#include "LBOBulkEnergy.h"
#include "TotalEnergyAll.h"

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
  registerAux(PolarizationVortexAux);
  registerAux(TensorPressureAux);
  registerAux(BulkEnergyDensity);
  registerAux(WallEnergyDensity);
  registerAux(ElecFieldAux);
  registerAux(ChernSimonsDensity);
  registerAux(ChernSimonsDensityMag);
  registerAux(CurlP);
  registerAux(ZCompCurlP);
  registerAux(CurlPMag);
  registerAux(BandGapAuxZnO);
  registerAux(BandGapAuxTiO2);
  registerAux(SurfaceChargeAux);
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
  registerAux(SemiconductingChargeCarriersAux);
  registerAux(Birefringence);
  registerAux(ChangeInRefractiveIndex);
  registerAux(RefractiveIndex);
  registerAux(ChangeInRefractiveIndexWithPolar);
  registerAux(ChangeInRefractiveIndexWithGCoeffPolar);
  registerAux(PkNorm);
  registerAux(SkyrmionChargeDensityZ);
  registerAux(SemiconductingChargeCarriersPolyLogAux);
  registerAux(ChangeInRefractiveIndexElectro);

  ///Kernels
  registerKernel(ModifiedStressDivergenceTensors);
  registerKernel(BulkEnergyDerivativeSixth);
  registerKernel(BulkEnergyDerivativeSixthAlt);
  registerKernel(BulkEnergyDerivativePSTO);
  registerKernel(RotatedBulkEnergyDerivativeSixth);
  registerKernel(RotatedBulkEnergyDerivative);
  registerKernel(NoStdBulkEnergyDerivativeSixth);
  registerKernel(BulkEnergyDerivativeSixthCoupledT);
  registerKernel(WallEnergyDerivative);
  registerKernel(WallEnergyDerivativeAlt);
  registerKernel(RotatedWallEnergyDerivative);
  registerKernel(TimeDerivativeScaled);
  registerKernel(FerroelectricCouplingP);
  registerKernel(FluctuationKernel);
  registerKernel(KarmanenkoDriver);
  registerKernel(KappaTDiffusion);
  registerKernel(AnisotropyEnergy);
  registerKernel(LBOBulkEnergyDeriv);
  registerKernel(DepolEnergy);
  registerKernel(ThomasFermiPotential);
  registerKernel(ThomasFermiTerm);
  registerKernel(AnisotropicElectrostatics);
  registerKernel(RenormalizedFreeEnergy);
  registerKernel(FerroelectricCouplingX);
  registerKernel(StressDivergenceTensorsScaled);
  registerKernel(PolarElectricEStrong);
  registerKernel(PolarElectricPStrong);
  registerKernel(Electrostatics);
  registerKernel(CoeffParamDiffusion);
  registerKernel(SemiconductorChargeCarriers);
  registerKernel(NerstPlanckDrivingTerm);
  registerKernel(NerstPlanckDiffusive);
  registerKernel(FreeChargeContribution);
  registerKernel(HoleChargeContribution);
  registerKernel(AcceptorIonContribution);
  registerKernel(SemiconductorChargeCarriersPolyLog);
  registerKernel(PolarElectricEStrongAlt);
  registerKernel(ConversePiezoelectricStrain);
  registerKernel(PiezoelectricStrainCharge);

  ///registerInterfaceKernels
  registerInterfaceKernel(InterfaceDiffusion);

  ///Postprocessors
  registerPostprocessor(BulkEnergy);
  registerPostprocessor(BulkEnergyPSTO);
  registerPostprocessor(BulkEnergyCoupledT);
  registerPostprocessor(WallEnergy);
  ///registerPostprocessor(ChernSimonsNumber);
  registerPostprocessor(ElectrostaticEnergy);
  registerPostprocessor(TotalEnergy);
  registerPostprocessor(TotalEnergyPSTO);
  registerPostprocessor(TotalEnergyPSTOcoupled);
  registerPostprocessor(CoupledEnergyPSTO);
  registerPostprocessor(TotalEnergyFlow);
  registerPostprocessor(TotalEnergyFlowNoElast);
  registerPostprocessor(TotalEnergyFlowNoElastNoElec);
  registerPostprocessor(ElasticEnergy);
  registerPostprocessor(CoupledEnergy);
  registerPostprocessor(ElectrostrictiveEnergy);
  registerPostprocessor(ThermalEnergy);
  registerPostprocessor(CoupledEnergyCheckShear);
  registerPostprocessor(GrainSize);
  registerPostprocessor(DepolarizationEnergy);
  registerPostprocessor(AnisotropicEnergy);
  registerPostprocessor(TotalEnergyG);
  registerPostprocessor(RenormalizedBulkEnergy);
  registerPostprocessor(TotalEnergyP);
  registerPostprocessor(TotalEnergySkFlow);
  registerPostprocessor(TotalWinding);
  registerPostprocessor(EnergyRatePostprocessor);
  registerPostprocessor(LBOBulkEnergy);
  registerPostprocessor(TotalEnergyAll);
  registerPostprocessor(ExtElectrostaticEnergy);

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

  ///InitialConditions
  registerInitialCondition(FluctuationsIC);

#undef registerObject
#define registerObject(name) factory.regLegacy<name>(stringifyName(name))
}
