/**
 * @file   EulerSkyrmionThetaTerm.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Feb. 20. 2017
 *
 */

#include "EulerSkyrmionThetaTerm.h"

class EulerSkyrmionThetaTerm;

template<>
InputParameters validParams<EulerSkyrmionThetaTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("theta", "The theta variable");
  params.addRequiredCoupledVar("P", "The polar magnitude variable");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("xi0", 1.0, "the domain wall coefficient");
  return params;
}

EulerSkyrmionThetaTerm::EulerSkyrmionThetaTerm(const InputParameters & parameters)
  :Kernel(parameters),
   _theta_var(coupled("theta")),
   _P_var(coupled("P")),
   _theta(coupledValue("theta")),
   _second_u(second()),
   _second_test(secondTest()),
   _second_phi(secondPhi()),
   _P(coupledValue("P")),
   _len_scale(getParam<Real>("len_scale")),
   _xi0(getParam<Real>("xi0"))
{
}

Real
EulerSkyrmionThetaTerm::computeQpResidual()
{
  return - _test[_i][_qp] * _xi0 * _xi0 * _P[_qp] * _second_u[_qp](0,0);
}

Real
EulerSkyrmionThetaTerm::computeQpJacobian()
{
  return - _test[_i][_qp] * _xi0 * _xi0 * _P[_qp] * _second_phi[_j][_qp](0,0);
}

Real
EulerSkyrmionThetaTerm::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _P_var)
    return - _test[_i][_qp] * _xi0 * _xi0 * _phi[_j][_qp] * _second_u[_qp](0,0);
  else
    return 0.0;
}
