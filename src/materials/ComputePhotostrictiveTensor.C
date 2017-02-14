/****************************************************************/
/* Computes a rank 4 electrostrictive tensor                    */
/****************************************************************/

#include "ComputeElasticityTensor.h"
#include "ComputeRotatedElasticityTensorBase.h"
#include "ComputePhotostrictiveTensor.h"
#include "RotationTensor.h"
#include "RankFourTensor.h"

template<>
InputParameters validParams<ComputePhotostrictiveTensor>()
{
  InputParameters params = validParams<ComputeRotatedPhotostrictiveTensorBase>();
  params.addClassDescription("Compute a photostrictive tensor.");
  params.addRequiredParam<std::vector<Real> >("P_mnkl", "elasto-optic tensor for material");
  params.addParam<MooseEnum>("fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  return params;
}

ComputePhotostrictiveTensor::ComputePhotostrictiveTensor(const InputParameters & parameters) :
    ComputeRotatedPhotostrictiveTensorBase(parameters),
    _Pmnkl(getParam<std::vector<Real> >("P_mnkl"), (RankFourTensor::FillMethod)(int)getParam<MooseEnum>("fill_method"))
{
  /// Define a rotation according to Euler angle parameters
  RotationTensor R(_Euler_angles); // R type: RealTensorValue
  /// rotate electrostrictive tensor -- note that it needs to be collinear with the elasticity tensor _always_
  _Pmnkl.rotate(R);
}

void
ComputePhotostrictiveTensor::computeQpPhotostrictiveTensor()
{
  ///Assign a photostrictive tensor at a given quad point. This will be reworked eventually for constant _qp.
  _photostrictive_tensor[_qp] = _Pmnkl;
}
