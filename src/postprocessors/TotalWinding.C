/**
 * @file   TotalWinding.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   
 *
 * @brief Integral over the volume of the winding number density q.
 *
 */

#include "TotalWinding.h"

template<>
InputParameters validParams<TotalWinding>()
{

  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("q", "The winding number density");
  return params;
}

TotalWinding::TotalWinding(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _q(coupledValue("q"))
{
}

Real
TotalWinding::computeQpIntegral()
{
  return _q[_qp];
}
