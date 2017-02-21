/**
 * @file   EulerSkyrmionPTempTerm.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Feb. 20. 2017
 *
 */

#include "EulerSkyrmionPTempTerm.h"

class EulerSkyrmionPTempTerm;

template<>
InputParameters validParams<EulerSkyrmionPTempTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("theta", "The theta variable");
  params.addRequiredCoupledVar("P", "The polar magnitude variable");
  params.addParam<Real>("t", 1.0, "the alpha1 temperature constant");
  params.addParam<Real>("kappa", 1.0, "the anisotropy constant");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("xi0", 1.0, "the domain wall coefficient");
  return params;
}

EulerSkyrmionPTempTerm::EulerSkyrmionPTempTerm(const InputParameters & parameters)
  :Kernel(parameters),
   _theta_var(coupled("theta")),
   _P_var(coupled("P")),
   _theta(coupledValue("theta")),
   _theta_grad(coupledGradient("theta")),
   _P(coupledValue("P")),
   _t(getParam<Real>("t")),
   _kappa(getParam<Real>("kappa")),
   _len_scale(getParam<Real>("len_scale")),
   _xi0(getParam<Real>("xi0"))
{
}

Real
EulerSkyrmionPTempTerm::computeQpResidual()
{
  return _test[_i][_qp] * (
  _t + _xi0 * _xi0 * _theta_grad[_qp](0) * _theta_grad[_qp](0)
  + (_kappa + _xi0 * _xi0 / (_q_point[_qp](0) * _q_point[_qp](0))) * std::sin(_theta[_qp]) * std::sin(_theta[_qp]) ) * _P[_qp];
}

Real
EulerSkyrmionPTempTerm::computeQpJacobian()
{
  return _test[_i][_qp] * (
  _t + _xi0 * _xi0 * _theta_grad[_qp](0) * _theta_grad[_qp](0)
  + (_kappa + _xi0 * _xi0 / (_q_point[_qp](0) * _q_point[_qp](0))) * std::sin(_theta[_qp]) * std::sin(_theta[_qp]) ) * _phi[_j][_qp];
}

Real
EulerSkyrmionPTempTerm::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
