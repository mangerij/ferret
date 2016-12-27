/**
 * @file   CubicDielectricTensor.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief calculate the components of the anisotropic dielectric tensor 
 *        assuming cubic symmetry of the parent phase
 *
 */
#include "CubicDielectricTensor.h"

template<>
InputParameters validParams<CubicDielectricTensor>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("alpha1", "alpha1 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha11", "alpha11 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha12", "alpha12 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha111", "alpha111 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha112", "alpha112 coefficient of the Landau expansion");
  params.addRequiredParam<Real>("alpha123", "alpha123 coefficient of the Landau expansion");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

CubicDielectricTensor::CubicDielectricTensor(const InputParameters & parameters) :
  AuxKernel(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _alpha1(getParam<Real>("alpha1")),
  _alpha11(getParam<Real>("alpha11")),
  _alpha12(getParam<Real>("alpha12")),
  _alpha111(getParam<Real>("alpha111")),
  _alpha112(getParam<Real>("alpha112")),
  _alpha123(getParam<Real>("alpha123")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
CubicDielectricTensor::computeValue()
{
  /**if (_component = 0 && _off_component == 0)
  {
    return 2 (\[Alpha]1 + 15 Px^4 \[Alpha]111 + Py^4 \[Alpha]112 + 
     Pz^4 \[Alpha]112 + 
     6 Px^2 (\[Alpha]11 + (Py^2 + Pz^2) \[Alpha]112) + Pz^2 \[Alpha]12 +
      Py^2 (\[Alpha]12 + Pz^2 \[Alpha]123) + (c11 q11 + 
        2 c12 q12) \[CurlyEpsilon][1, 
       1] + (c12 q11 + c11 q12 + c12 q12) \[CurlyEpsilon][2, 
       2] + (c12 q11 + c11 q12 + c12 q12) \[CurlyEpsilon][3, 3]);
  }
  else if ()
  {}

  else if ()
  {}

  else if ()
  {}

  else if ()
  {}
  else if ()
  {}

  else ()
  {} **/
  return 0.0;
}
