#include "BoundCharge.h"

template<>

InputParameters validParams<BoundCharge>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
  params.addRequiredCoupledVar("potential_ext", "The external electric potential variable");
  return params;
}


BoundCharge::BoundCharge(const InputParameters & parameters) :
  AuxKernel(parameters),
   _potential_int_grad(coupledGradient("potential_int")),
   _potential_ext_grad(coupledGradient("potential_ext"))
{
}

Real
BoundCharge::computeValue()

{
    return _potential_ext_grad[_qp](0)+_potential_ext_grad[_qp](1)+_potential_ext_grad[_qp](2)+_potential_int_grad[_qp](0)+_potential_int_grad[_qp](1)+_potential_int_grad[_qp](2);
}
