#include "PolarizationComponentValue.h"

template<>
InputParameters validParams<PolarizationComponentValue>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar", "Component of the polarization");

  return params;
}

PolarizationComponentValue::PolarizationComponentValue(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _polar(coupledValue("polar"))
{
}

Real
PolarizationComponentValue::computeQpIntegral()
{
  return _polar[_qp] * _polar[_qp];
}
