/**
 * @file   WindingNumberDensity.C
 *
 * @brief calculate the winding number density in 3D for the polarization field P
 *
 */
#include "WindingNumberDensity.h"

template<>
InputParameters validParams<WindingNumberDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return params;
}

WindingNumberDensity::WindingNumberDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z"))
{
}

Real
WindingNumberDensity::computeValue()
{
  return
  _polar_x[_qp] * (_polar_y_grad[_qp](0) * _polar_x_grad[_qp](1) * _polar_y_grad[_qp](2) - _polar_x_grad[_qp](0) *_polar_y_grad[_qp](1) *_polar_y_grad[_qp](2) + _polar_z_grad[_qp](0) *_polar_x_grad[_qp](1) *_polar_z_grad[_qp](2) - _polar_x_grad[_qp](0) *_polar_z_grad[_qp](1) *_polar_z_grad[_qp](2))
  + _polar_y[_qp] * (-_polar_y_grad[_qp](0) * _polar_x_grad[_qp](1) * _polar_x_grad[_qp](2) + _polar_x_grad[_qp](0) *_polar_y_grad[_qp](1) *_polar_x_grad[_qp](2) + _polar_z_grad[_qp](0) *_polar_y_grad[_qp](1) *_polar_z_grad[_qp](2) - _polar_y_grad[_qp](0) *_polar_z_grad[_qp](1) *_polar_z_grad[_qp](2))
  + _polar_z[_qp] * (-_polar_z_grad[_qp](0) * _polar_x_grad[_qp](1) * _polar_x_grad[_qp](2) + _polar_x_grad[_qp](0) *_polar_z_grad[_qp](1) *_polar_x_grad[_qp](2) - _polar_z_grad[_qp](0) *_polar_y_grad[_qp](1) *_polar_y_grad[_qp](2) + _polar_y_grad[_qp](0) *_polar_z_grad[_qp](1) *_polar_y_grad[_qp](2));
}
