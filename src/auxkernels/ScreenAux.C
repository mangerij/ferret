#include "ScreenAux.h"

template<>

InputParameters validParams<ScreenAux>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("potential_int", "The internal potential");
  return params;
}


ScreenAux::ScreenAux( const std::string & name, InputParameters parameters ) :
  AuxKernel( name, parameters ),
  _normals(_var.normals()),
  _potential_int_grad(coupledGradient("potential_int"))
{
}

Real
ScreenAux::computeValue()
{
  RealVectorValue n(_normals[_qp]);
  return _potential_int_grad[_qp](0) * n(0) + _potential_int_grad[_qp](1) * n(1) + _potential_int_grad[_qp](2) * n(2);
}
