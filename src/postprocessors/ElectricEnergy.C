/**
 * @file   ElectricEnergy.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jul 30 17:09:44 2013
 *
 * @brief
 *
 *
 */

#include "ElectricEnergy.h"

template<>
InputParameters validParams<ElectricEnergy>()
{
  //TODO: inherit from an appropriate postprocessor
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addRequiredCoupledVar("potential", "The electric potential");
  params.addRequiredParam<Real>("permittivity", "permittivity");
  return params;
}

ElectricEnergy::ElectricEnergy(const std::string & name, InputParameters parameters) :
  ElementIntegralPostprocessor(name, parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _potential_grad(coupledGradient("potential")),
  _permittivity(getParam<Real>("permittivity"))
{
}

Real
ElectricEnergy::computeQpIntegral()
{
  RealVectorValue E, D, P;
  P(0)=_polar_x[_qp];P(1)=_polar_y[_qp];P(2)=_polar_z[_qp];
  E=_potential_grad[_qp]*(-1.0);
  D=E*_permittivity-P;
  return 0.5*D*E;
}
