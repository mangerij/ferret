/* This class is solely for testing electrostatics. We implement just one potential
with two different dielectric permittivities on each block and then take the solution potential
with the definition of the electric displacement to calculate P */

#include "PzFieldAux.h"

template<>

InputParameters validParams<PzFieldAux>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<Real>("permittivity_int", "internal permittivity");
  params.addRequiredParam<Real>("permittivity_ext", "external permittivity");
//  params.addRequiredCoupledVar("potential_int", "The internal electric potential variable");
  params.addRequiredCoupledVar("potential_ext", "The external electric potential variable");
  return params;
}


PzFieldAux::PzFieldAux( const std::string & name, InputParameters parameters ) :
  AuxKernel( name, parameters ),
   _permittivity_int(getParam<Real>("permittivity_int")),
   _permittivity_ext(getParam<Real>("permittivity_ext")),
//   _potential_int_grad(coupledGradient("potential_int")),
   _potential_ext_grad(coupledGradient("potential_ext"))
{
}

Real
PzFieldAux::computeValue()

{
  return - ( _permittivity_int - _permittivity_ext) * _potential_ext_grad[_qp](2);
}
