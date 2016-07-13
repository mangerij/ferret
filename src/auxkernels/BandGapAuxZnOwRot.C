#include "BandGapAuxZnOwRot.h"
#include "TensorMechanicsMaterial.h"
#include "RotationTensor.h"
#include "Material.h"
#include "RankTwoTensor.h"

template<>

InputParameters validParams<BandGapAuxZnOwRot>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("relaxed_energy", 0.0,"relaxed energy");
  params.addParam<Real>("biaxial_strain_rate", 0.0, "uniaxial strain rate");
  params.addParam<Real>("uniaxial_strain_rate", 0.0, "biaxial strain rate");
  params.addParam<Real>("biaxial_relaxation_coeff", 0.0, "biaxial relaxation coeff");
  params.addParam<Real>("poisson_ratio", 0.0, "Poisson ratio");
  params.addParam<Real>("euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("euler_angle_3", 0.0, "Euler angle in direction 3");
  return params;
}


BandGapAuxZnOwRot::BandGapAuxZnOwRot(const InputParameters & parameters) :
    AuxKernel(parameters),
    _strain(getMaterialProperty<RankTwoTensor>("elastic_strain")),
    _du(getParam<Real>("biaxial_strain_rate")),
    _db(getParam<Real>("uniaxial_strain_rate")),
    _E0(getParam<Real>("relaxed_energy")),
    _Rb(getParam<Real>("biaxial_relaxation_coeff")),
    _nu(getParam<Real>("poisson_ratio")),
    _Euler_angles(getParam<Real>("euler_angle_1"),
                  getParam<Real>("euler_angle_2"),
                  getParam<Real>("euler_angle_3"))
{
  RotationTensor R(_Euler_angles);
}

Real
BandGapAuxZnOwRot::computeValue()

{
  RotationTensor R(_Euler_angles);
    return _E0 + (1/(1-_Rb))*((_db+_du*_Rb)*0.5*(((R(0,0)*_strain[_qp](0,0)*R(0,0)+R(0,1)*_strain[_qp](1,0)*R(0,0)+R(0,2)*_strain[_qp](2,0)*R(0,0)+
    R(0,0)*_strain[_qp](0,2)*R(2,0)+R(0,1)*_strain[_qp](1,2)*R(2,0)+R(0,2)*_strain[_qp](2,2)*R(2,0)))+
    R(1,0)*_strain[_qp](0,0)*R(0,1)+R(1,1)*_strain[_qp](1,0)*R(0,1)+R(1,2)*_strain[_qp](2,0)*R(0,1)+
    R(1,0)*_strain[_qp](0,1)*R(1,1)+R(1,1)*_strain[_qp](1,1)*R(1,1)+R(1,2)*_strain[_qp](2,1)*R(1,1)+
    R(1,0)*_strain[_qp](0,2)*R(2,1)+R(1,1)*_strain[_qp](1,2)*R(2,1)+R(1,2)*_strain[_qp](2,2)*R(2,1))
    +(_du+_nu*_db)*(R(2,0)*_strain[_qp](0,0)*R(0,2)+R(2,1)*_strain[_qp](1,0)*R(0,2)+R(2,2)*_strain[_qp](2,0)*R(0,2)+
    R(2,0)*_strain[_qp](0,1)*R(1,2)+R(2,1)*_strain[_qp](1,1)*R(1,2)+R(2,2)*_strain[_qp](2,1)*R(1,2)+
    R(2,0)*_strain[_qp](0,2)*R(2,2)+R(2,1)*_strain[_qp](1,2)*R(2,2)+R(2,2)*_strain[_qp](2,2)*R(2,2)));
}
