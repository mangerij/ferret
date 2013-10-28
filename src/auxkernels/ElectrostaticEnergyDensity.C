/**
 * @file   ElectrostaticEnergyDensity.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Oct  2 17:14:17 2013
 *
 * @brief
 *
 *
 */


#include "ElectrostaticEnergyDensity.h"

template<>
InputParameters validParams<ElectrostaticEnergyDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("potential_int", "The internal electrostatic potential");
  params.addRequiredCoupledVar("potential_ext", "The external electrostatic potential");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

ElectrostaticEnergyDensity::ElectrostaticEnergyDensity(const std::string & name, InputParameters parameters) :
  AuxKernel(name, parameters),
  _potential_int_grad(coupledGradient("potential_int")),
  _potential_ext_grad(coupledGradient("potential_ext")),
  _polar_x(coupledValueOld("polar_x")),
  _polar_y(coupledValueOld("polar_y")),
  _polar_z(coupledValueOld("polar_z")),
  _len_scale(getParam<Real>("len_scale"))
{}


Real
ElectrostaticEnergyDensity::computeValue()
{
  RealVectorValue P;
  P(0)=_polar_x[_qp];P(1)=_polar_y[_qp];P(2)=_polar_z[_qp];
  return 0.5*(P*_potential_int_grad[_qp])*pow(_len_scale,2.0)+(P*_potential_ext_grad[_qp])*pow(_len_scale,2.0);
}
