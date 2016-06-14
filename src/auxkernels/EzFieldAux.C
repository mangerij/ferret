#include "EzFieldAux.h"

template<>

InputParameters validParams<EzFieldAux>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addCoupledVar("potential_int", "The internal electric potential variable");
  params.addCoupledVar("potential_ext", "The external electric potential variable");
  return params;
}


EzFieldAux::EzFieldAux(const InputParameters & parameters) :
  AuxKernel(parameters ),
   _potential_int_grad(coupledGradient("potential_int")),
   _potential_ext_grad(coupledGradient("potential_ext"))
{
}

Real
EzFieldAux::computeValue()

{
    return - _potential_int_grad[_qp](2) - _potential_ext_grad[_qp](2);
}
