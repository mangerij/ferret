/**
 * @file   ThermalEnergy.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 */

#include "ThermalEnergy.h"

template<>
InputParameters validParams<ThermalEnergy>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("temperature", "The local temperature");
  params.addParam<Real>("const", 1.0, "thermal relation constant (e.g. Boltzmann constant)");
  return params;
}

ThermalEnergy::ThermalEnergy(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _temperature(coupledValue("temperature")),
   _const(getParam<Real>("const"))
{
}

Real
ThermalEnergy::computeQpIntegral()
{
  return (3.0/2.0) * _const * _temperature[_qp];
}
