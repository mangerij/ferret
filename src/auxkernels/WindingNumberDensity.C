/**
 * @file   WindingNumberDensity.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief calculate the winding number density in 3D for the polarization field P with 
 *        out-of-plane direction along z (this is the "t" coordinate system in 
 *        V. Stepkova and J. Hlinka, Phase Transitions, 2017 Vol. 90, No. 1, 11-16)
 */

#include "WindingNumberDensity.h"

template<>
InputParameters validParams<WindingNumberDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("norm_polar_x", "The x component of the normalized polarization");
  params.addRequiredCoupledVar("norm_polar_y", "The y component of the normalized polarization");
  params.addRequiredCoupledVar("norm_polar_z", "The z component of the normalized polarization");
  return params;
}

WindingNumberDensity::WindingNumberDensity(const InputParameters & parameters) :
  AuxKernel(parameters),
  _norm_polar_x(coupledValue("norm_polar_x")),
  _norm_polar_y(coupledValue("norm_polar_y")),
  _norm_polar_z(coupledValue("norm_polar_z")),
  _norm_polar_x_grad(coupledGradient("norm_polar_x")),
  _norm_polar_y_grad(coupledGradient("norm_polar_y")),
  _norm_polar_z_grad(coupledGradient("norm_polar_z"))
{
}

Real
WindingNumberDensity::computeValue()
{
  // Pz (-xPy yPx + xPx yPy) + Py (xPz yPx - xPx yPz) +  Px (-xPz yPy + xPy yPz)
  return (1.0 / (4.0 * 3.14159265359)) * (
   _norm_polar_x[_qp] * (-_norm_polar_y_grad[_qp](0) * _norm_polar_x_grad[_qp](1) + _norm_polar_x_grad[_qp](0) * _norm_polar_y_grad[_qp](1))
 + _norm_polar_y[_qp] * (_norm_polar_z_grad[_qp](0) * _norm_polar_x_grad[_qp](1) - _norm_polar_x_grad[_qp](0) * _norm_polar_z_grad[_qp](1))
 + _norm_polar_z[_qp] * (-_norm_polar_z_grad[_qp](0) * _norm_polar_y_grad[_qp](1) + _norm_polar_y_grad[_qp](0) * _norm_polar_z_grad[_qp](1)));
}
