#include "ExFieldAux.h"

template<>

InputParameters validParams<ExFieldAux>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
  params.addRequiredCoupledVar("potential_ext", "The external electric potential variable");
  return params;
}


ExFieldAux::ExFieldAux(const InputParameters & parameters) :
  AuxKernel(parameters),
   _potential_int_grad(coupledGradient("potential_int")),
   _potential_ext_grad(coupledGradient("potential_ext"))
{
}

Real
ExFieldAux::computeValue()

{
    return - _potential_int_grad[_qp](0) - _potential_ext_grad[_qp](0);
}
