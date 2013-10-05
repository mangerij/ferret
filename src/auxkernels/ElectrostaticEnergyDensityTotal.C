/**
 * @file   ElectrostaticEnergyDensityTotal.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Oct  2 17:30:58 2013
 */
#include "ElectrostaticEnergyDensityTotal.h"

template<>
InputParameters validParams<ElectrostaticEnergyDensityTotal>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<Real>("permittivity", "permittivity");
  params.addRequiredCoupledVar("potential", "The electrostatic potential");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return params;
}

ElectrostaticEnergyDensityTotal::ElectrostaticEnergyDensityTotal(const std::string & name, InputParameters parameters) :
  AuxKernel(name, parameters),
  _permittivity(getParam<Real>("permittivity")),
  _potential_grad(coupledGradient("potential")),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z"))
{}

Real
ElectrostaticEnergyDensityTotal::computeValue()
{
  return _permittivity*_potential_grad[_qp].size_sq()+_potential_grad[_qp](0)*_polar_x[_qp]+ _potential_grad[_qp](1)*_polar_y[_qp]+ _potential_grad[_qp](2)*_polar_z[_qp];
}
