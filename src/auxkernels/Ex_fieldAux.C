#include "Ex_fieldAux.h"

template<>

InputParameters validParams<Ex_fieldAux>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
  params.addRequiredCoupledVar("potential_ext", "The external electric potential variable");
  return params;
}


Ex_fieldAux::Ex_fieldAux( const std::string & name, InputParameters parameters ) :
  AuxKernel( name, parameters ),
   _potential_int_grad(coupledGradient("potential_int")),
   _potential_ext_grad(coupledGradient("potential_ext"))
{
}

Real
Ex_fieldAux::computeValue()

{
    return -_potential_int_grad[_qp](0)-_potential_ext_grad[_qp](0);
}


