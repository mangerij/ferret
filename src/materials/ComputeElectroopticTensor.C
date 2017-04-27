/****************************************************************/
/* Computes a rank 4 electrostrictive tensor                    */
/****************************************************************/


#include "ComputeElectroopticTensor.h"
#include "RotationTensor.h"
#include "RankThreeTensor.h"

template<>
InputParameters validParams<ComputeElectroopticTensor>()
{
  InputParameters params = validParams<ComputeRotatedElectroopticTensorBase>();
  params.addClassDescription("Compute an electrooptic tensor.");
  params.addRequiredParam<std::vector<Real> >("r_ijk", "electrooptic tensor for material");
  params.addParam<MooseEnum>("fill_method", RankThreeTensor::fillMethodEnum() = "general", "The fill method");
  return params;
}

ComputeElectroopticTensor::ComputeElectroopticTensor(const InputParameters & parameters) :
    ComputeRotatedElectroopticTensorBase(parameters),
    _rijk(getParam<std::vector<Real> >("r_ijk"), (RankThreeTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue TODO: DOUBLE CHECK THAT THIS INDEED DOES WORK.
  /// rotate electrostrictive tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _rijk.rotate(R);
}

void
ComputeElectroopticTensor::computeQpElectroopticTensor()
{
  ///Assign a photostrictive tensor at a given quad point. This will be reworked eventually for constant _qp.
  _electrooptic_tensor[_qp] = _rijk;
}

