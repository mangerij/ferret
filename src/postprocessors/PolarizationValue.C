#include "PolarizationValue.h"
#include <cmath>

template<>
InputParameters validParams<PolarizationValue>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");

  return params;
}

PolarizationValue::PolarizationValue(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z"))
{
}

Real
PolarizationValue::computeQpIntegral()
{
  return std::sqrt( (_polar_x[_qp] * _polar_x[_qp]) + (_polar_y[_qp] * _polar_y[_qp]) + (_polar_z[_qp] * _polar_z[_qp]) );
}
