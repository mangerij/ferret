#include "BandGapAuxTiO2.h"

template<>

InputParameters validParams<BandGapAuxTiO2>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("relaxed_energy", 0.0,"relaxed energy");
  params.addParam<Real>("biaxial_stress_rate", 0.0, "biaxial stress rate");
  params.addParam<Real>("uniaxial_stress_rate", 0.0, "uniaxial stress rate");
  return params;
}


BandGapAuxTiO2::BandGapAuxTiO2(const InputParameters & parameters) :
  AuxKernel(parameters),
   _stress(getMaterialProperty<RankTwoTensor>("stress")),
   _ba(getParam<Real>("biaxial_stress_rate")),
   _bc(getParam<Real>("uniaxial_stress_rate")),
   _E0(getParam<Real>("relaxed_energy"))
{
}

Real
BandGapAuxTiO2::computeValue()

{
    return _E0 + _ba*_stress[_qp](0,0)+_bc*_stress[_qp](2,2);
}


