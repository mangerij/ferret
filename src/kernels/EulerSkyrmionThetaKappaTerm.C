/**
 * @file   EulerSkyrmionThetaKappaTerm.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Feb. 20. 2017
 *
 */

#include "EulerSkyrmionThetaKappaTerm.h"

class EulerSkyrmionThetaKappaTerm;

template<>
InputParameters validParams<EulerSkyrmionThetaKappaTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("theta", "The theta variable");
  params.addRequiredCoupledVar("P", "The polar magnitude variable");
  params.addParam<Real>("kappa", 1.0, "the anisotropy constant kappa");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("xi0", 1.0, "the domain wall coefficient");
  return params;
}

EulerSkyrmionThetaKappaTerm::EulerSkyrmionThetaKappaTerm(const InputParameters & parameters)
  :Kernel(parameters),
   _theta_var(coupled("theta")),
   _P_var(coupled("P")),
   _theta(coupledValue("theta")),
   _second_u(second()),
   _second_test(secondTest()),
   _second_phi(secondPhi()),
   _P(coupledValue("P")),
   _kappa(getParam<Real>("kappa")),
   _len_scale(getParam<Real>("len_scale")),
   _xi0(getParam<Real>("xi0"))
{
}

Real
EulerSkyrmionThetaKappaTerm::computeQpResidual()
{
  return _test[_i][_qp] * 0.5 * (_kappa + _xi0 * _xi0 / (_q_point[_qp](0) * _q_point[_qp](0))) * _P[_qp] * std::sin(2.0 * _theta[_qp]);
}

Real
EulerSkyrmionThetaKappaTerm::computeQpJacobian()
{
  return _test[_i][_qp] * _phi[_j][_qp] * (_kappa + _xi0 * _xi0 / (_q_point[_qp](0) * _q_point[_qp](0))) * _P[_qp] * std::cos(2.0 * _theta[_qp]);
}

Real
EulerSkyrmionThetaKappaTerm::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _P_var)
    return _test[_i][_qp] * _phi[_j][_qp] * 0.5 * (_kappa + _xi0 * _xi0 / (_q_point[_qp](0) * _q_point[_qp](0)))  * std::sin(2.0 * _theta[_qp]);
  else
    return 0.0;
}
