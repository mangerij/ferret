/**
 * @file   ComputeRotatedElectroopticTensorBase.C
 * @author J. Mangeri <john.mangeri@uconn.edu
 *
 * @brief  Base class for linear electrooptic material
 *         Note that the electrooptic tensor r_{ijk}
 *         must be collinear with C_{ijkl}
 */

#include "ComputeRotatedElectroopticTensorBase.h"
#include "RotationTensor.h"

template<>
InputParameters validParams<ComputeRotatedElectroopticTensorBase>()
{
  InputParameters params = validParams<ComputeElectroopticTensorBase>();
  params.addParam<Real>("euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("euler_angle_3", 0.0, "Euler angle in direction 3");
  return params;
}

ComputeRotatedElectroopticTensorBase::ComputeRotatedElectroopticTensorBase(const InputParameters & parameters) :
    ComputeElectroopticTensorBase(parameters),
    _Euler_angles(getParam<Real>("euler_angle_1"),
                  getParam<Real>("euler_angle_2"),
                  getParam<Real>("euler_angle_3"))
{
}
