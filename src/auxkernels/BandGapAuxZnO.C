#include "BandGapAuxZnO.h"
#include "TensorMechanicsMaterial.h"

template<>

InputParameters validParams<BandGapAuxZnO>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("relaxed_energy", 0.0,"relaxed energy");
  params.addParam<Real>("biaxial_strain_rate", 0.0, "uniaxial strain rate");
  params.addParam<Real>("uniaxial_strain_rate", 0.0, "biaxial strain rate");
  params.addParam<Real>("biaxial_relaxation_coeff", 0.0, "biaxial relaxation coeff");
  params.addParam<Real>("poisson_ratio", 0.0, "Poisson ratio");
  return params;
}


BandGapAuxZnO::BandGapAuxZnO(const InputParameters & parameters) :
  AuxKernel(parameters),
   _strain(getMaterialProperty<RankTwoTensor>("elastic_strain")),
   _du(getParam<Real>("biaxial_strain_rate")),
   _db(getParam<Real>("uniaxial_strain_rate")),
   _E0(getParam<Real>("relaxed_energy")),
   _Rb(getParam<Real>("biaxial_relaxation_coeff")),
   _nu(getParam<Real>("poisson_ratio"))
{
}

Real
BandGapAuxZnO::computeValue()

{
    return _E0 + (1/(1-_Rb))*((_db+_du*_Rb)*0.5*(_strain[_qp](1,1)+_strain[_qp](0,0))+(_du+_nu*_db)*_strain[_qp](2,2));
}


