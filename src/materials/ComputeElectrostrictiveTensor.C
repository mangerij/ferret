/****************************************************************/
/* Computes a rank 4 electrostrictive tensor
/*
/****************************************************************/

#include "ComputeElasticityTensor.h"
#include "ComputeRotatedElasticityTensorBase.h"
#include "ComputeElectrostrictiveTensor.h"
#include "RotationTensor.h"
#include "RankFourTensor.h"

template<>
InputParameters validParams<ComputeElectrostrictiveTensor>()
{
  InputParameters params = validParams<ComputeRotatedElectrostrictiveTensorBase>();
  params.addClassDescription("Compute an electrostrictive tensor.");
  params.addRequiredParam<std::vector<Real> >("Q_mnkl", "electrostrictive tensor for material");
  params.addRequiredParam<std::vector<Real> >("C_ijkl", "elastic stiffness tensor for material");
  params.addParam<MooseEnum>("fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  return params;
}

ComputeElectrostrictiveTensor::ComputeElectrostrictiveTensor(const InputParameters & parameters) :
    ComputeRotatedElectrostrictiveTensorBase(parameters),
    _Qmnkl(getParam<std::vector<Real> >("Q_mnkl"), (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method")),
    _Cijkl(getParam<std::vector<Real> >("C_ijkl"), (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor"))
{
  // Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue
  // rotate electrostrictive tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _Qmnkl.rotate(R);
  //contractions using namespace
  _qijkl = ElectrostrictiveTensorTools::computeProduct(_Cijkl, _Qmnkl);
  _QQijkl = ElectrostrictiveTensorTools::computeProduct(_qijkl, _Qmnkl); //contract again
}

void
ComputeElectrostrictiveTensor::computeQpElectrostrictiveTensor()
{
  //Assign an electrostrictive tensor at a given quad point -- in principle we DON'T want this?
  _electrostrictive_tensor[_qp] = _qijkl;
  _electrostrictive_tensorQ[_qp] = _QQijkl;
}
