#include "ComputeRotatedElectricalConductivityTensorBase.h"
#include "RotationTensor.h"
template <>
InputParameters
validParams<ComputeRotatedElectricalConductivityTensorBase>()
{
  InputParameters params = validParams<ComputeElectricalConductivityTensorBase>();
  params.addParam<Real>("euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("euler_angle_3", 0.0, "Euler angle in direction 3");
  return params;
}
ComputeRotatedElectricalConductivityTensorBase::ComputeRotatedElectricalConductivityTensorBase(
    const InputParameters & parameters)
  : ComputeElectricalConductivityTensorBase(parameters),
    _Euler_angles(getParam<Real>("euler_angle_1"),
                  getParam<Real>("euler_angle_2"),
                  getParam<Real>("euler_angle_3"))
{
}
