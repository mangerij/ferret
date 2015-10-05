#include "EyFieldAux.h"

template<>

InputParameters validParams<EyFieldAux>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
  params.addRequiredCoupledVar("potential_ext", "The external electric potential variable");
  return params;
}


EyFieldAux::EyFieldAux(const InputParameters & parameters) :
  AuxKernel(parameters ),
   _potential_int_grad(coupledGradient("potential_int")),
   _potential_ext_grad(coupledGradient("potential_ext"))
{
}

Real
EyFieldAux::computeValue()

{
    return -_potential_int_grad[_qp](1)-_potential_ext_grad[_qp](1);
}
