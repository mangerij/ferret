/*****************************************************/
/* Simply rotates the photostrictive tensor Pmnkl. */
/* Must be colinear with Cijkl                       */
/*****************************************************/

#include "ComputeRotatedPhotostrictiveTensorBase.h"
#include "RotationTensor.h"

template<>
InputParameters validParams<ComputeRotatedPhotostrictiveTensorBase>()
{
  InputParameters params = validParams<ComputePhotostrictiveTensorBase>();
  params.addParam<Real>("euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("euler_angle_3", 0.0, "Euler angle in direction 3");
  return params;
}

ComputeRotatedPhotostrictiveTensorBase::ComputeRotatedPhotostrictiveTensorBase(const InputParameters & parameters) :
    ComputePhotostrictiveTensorBase(parameters),
    _Euler_angles(getParam<Real>("euler_angle_1"),
                  getParam<Real>("euler_angle_2"),
                  getParam<Real>("euler_angle_3"))
{
}
