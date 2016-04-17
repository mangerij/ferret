/****************************************************************/
/* Computes a rank 4 electrostrictive tensor
/*
/****************************************************************/

#include "ComputeElectrostrictiveTensor.h"
#include "RotationTensor.h"

template<>
InputParameters validParams<ComputeElectrostrictiveTensor>()
{
  InputParameters params = validParams<ComputeRotatedElectrostrictiveTensorBase>();
  params.addClassDescription("Compute an electrostrictive tensor.");
  params.addRequiredParam<std::vector<Real> >("Q_mnkl", "electrostrictive tensor for material");
  params.addParam<MooseEnum>("fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  return params;
}

ComputeElectrostrictiveTensor::ComputeElectrostrictiveTensor(const InputParameters & parameters) :
    ComputeRotatedElectrostrictiveTensorBase(parameters),
    _Qmnkl(getParam<std::vector<Real> >("Q_mnkl"), (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  // Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue
  // rotate electrostrictive tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _Qmnkl.rotate(R);
}

void
ComputeElectrostrictiveTensor::computeQpElectrostrictiveCoefficients()
{
  _electrostrictivecoefficients[_qp] = _Qmnkl;
}


void
ComputeElectrostrictiveTensor::computeQpElectrostrictiveTensor()

{
  //Assign an electrostrictive tensor at a given quad point
  _electrostrictive_tensor[_qp].computeProduct(_Cijkl,_Qmnkl);
}
