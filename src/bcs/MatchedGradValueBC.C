/*************************************************************************
*  attempt to implement thin film epitaxy
****************************************************************/

#include "MatchedGradValueBC.h"

template<>
InputParameters validParams<MatchedGradValueBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredCoupledVar("u", "The variable which we are constructing a residual contribution for");
  params.addRequiredCoupledVar("v", "The variable whose value we are to match.");
  params.addRequiredParam<int>("component","Which component(0 for x, 1 for y, 2 for z) in traction is used");
  return params;
}

MatchedGradValueBC::MatchedGradValueBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _u_grad(coupledGradient("u")),
    _v(coupledValue("v")),
    _v_num(coupled("v")),
    _component(getParam<int>("component"))
{
}

Real
MatchedGradValueBC::computeQpResidual()
{
  return (_u_grad[_qp](_component) - _v[_qp])*_test[_i][_qp];
}

