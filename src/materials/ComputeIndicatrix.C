/**
 * @file   ComputeIndicatrix.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate the indicatrix (should be collinear 
 * with the crystal
 *
 *
 */

#include "ComputeIndicatrix.h"
#include "RotationTensor.h"

template<>
InputParameters validParams<ComputeIndicatrix>()
{
  InputParameters params = validParams<ComputeIndicatrixBase>();
  params.addClassDescription("Compute the indicatrix.");
  params.addParam<Real>("refractive_index_bulk_ordinary", 1.0,"refractive index along the ordinary axis");
  params.addParam<Real>("refractive_index_bulk_extraordinary", 1.0,"refractive index along the extraordinary axis");
  params.addParam<Real>("euler_angle_1", 0.0, "Euler angle in direction 1");
  params.addParam<Real>("euler_angle_2", 0.0, "Euler angle in direction 2");
  params.addParam<Real>("euler_angle_3", 0.0, "Euler angle in direction 3");
  return params;
}

ComputeIndicatrix::ComputeIndicatrix(const InputParameters & parameters) :
   ComputeIndicatrixBase(parameters),
   _no(getParam<Real>("refractive_index_bulk_ordinary")),
   _ne(getParam<Real>("refractive_index_bulk_extraordinary")),
   _Euler_angles(getParam<Real>("euler_angle_1"),
                 getParam<Real>("euler_angle_2"),
                 getParam<Real>("euler_angle_3"))
{
}

void
ComputeIndicatrix::computeQpIndicatrix()
{
  RotationTensor R(_Euler_angles);
  RealVectorValue n(_no, _no, _ne);
  for (unsigned int i = 0; i < 3; ++i)
    _indicatrix[_qp](i) = R(i, 0) * n(0) + R(1, i) * n(1) + R(i, 2) * n(2);
}

