#include "ComputeRotatedThermalConductivityTensorBase.h"
#include "RotationTensor.h"
template <>
InputParameters
validParams<ComputeRotatedThermalConductivityTensorBase>()
{
  InputParameters params = validParams<ComputeThermalConductivityTensorBase>();
  params.addParam<Real>("euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("euler_angle_3", 0.0, "Euler angle in direction 3");
  return params;
}
ComputeRotatedThermalConductivityTensorBase::ComputeRotatedThermalConductivityTensorBase(
    const InputParameters & parameters)
  : ComputeThermalConductivityTensorBase(parameters),
    _Euler_angles(getParam<Real>("euler_angle_1"),
                  getParam<Real>("euler_angle_2"),
                  getParam<Real>("euler_angle_3"))
{
}
