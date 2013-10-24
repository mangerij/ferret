/**
 * @file   ElectrostaticEnergyDensityCross.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Oct  2 17:14:17 2013
 *
 * @brief
 *
 *
 */


#include "ElectrostaticEnergyDensityCross.h"

template<>
InputParameters validParams<ElectrostaticEnergyDensityCross>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("potential_int", "The internal electrostatic potential");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

ElectrostaticEnergyDensityCross::ElectrostaticEnergyDensityCross(const std::string & name, InputParameters parameters) :
  AuxKernel(name, parameters),
  _potential_grad(coupledGradient("potential_int")),
  _polar_x(coupledValueOld("polar_x")),
  _polar_y(coupledValueOld("polar_y")),
  _polar_z(coupledValueOld("polar_z")),
  _len_scale(getParam<Real>("len_scale"))
{}


Real
ElectrostaticEnergyDensityCross::computeValue()
{
  return -0.5*(_potential_grad[_qp](0)*_polar_x[_qp]+ _potential_grad[_qp](1)*_polar_y[_qp]+ _potential_grad[_qp](2)*_polar_z[_qp])*pow(_len_scale,2.0);
}
