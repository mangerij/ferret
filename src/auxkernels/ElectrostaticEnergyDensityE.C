/**
 * @file   ElectrostaticEnergyDensityE.C   ElectrostaticEnergyDensityE.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Oct  2 15:06:12 2013   Wed Oct  2 15:05:11 2013
 * @brief Calculate ElectrostaticEnergyDensity |\grad\phi|^2
 */

#include "ElectrostaticEnergyDensityE.h"

template<>
InputParameters validParams<ElectrostaticEnergyDensityE>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("potential", "The electrostatic potential");
  return params;
}

ElectrostaticEnergyDensityE::ElectrostaticEnergyDensityE(const std::string & name, InputParameters parameters) :
  AuxKernel(name, parameters),
  _potential_grad(coupledGradient("potential"))
{}

Real
ElectrostaticEnergyDensityE::computeValue()
{
  return _potential_grad[_qp].size_sq();
}
