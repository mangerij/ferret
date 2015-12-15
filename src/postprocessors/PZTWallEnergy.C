/**
 * @file   PZTWallEnergy.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Dec 15 2015
 *
 * @brief
 *
 *
 */

#include "PZTWallEnergy.h"

template<>
InputParameters validParams<PZTWallEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addRequiredParam<Real>("G110","Domain wall penalty coefficients");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

PZTWallEnergy::PZTWallEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _polar_x_grad(coupledGradient("polar_x")),
  _polar_y_grad(coupledGradient("polar_y")),
  _polar_z_grad(coupledGradient("polar_z")),
  _G110(getParam<Real>("G110")),
  _len_scale(getParam<Real>("len_scale"))
{}

Real
PZTWallEnergy::computeQpIntegral()
{
  return (0.5*_G110*(pow(_polar_x_grad[_qp](0),2)+pow(_polar_y_grad[_qp](1),2)+pow(_polar_z_grad[_qp](2),2)) )*_len_scale;
}
