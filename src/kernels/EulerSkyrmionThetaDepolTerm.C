/**
 * @file   EulerSkyrmionThetaDepolTerm.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Feb. 20. 2017
 *
 */

#include "EulerSkyrmionThetaDepolTerm.h"

class EulerSkyrmionThetaDepolTerm;

template<>
InputParameters validParams<EulerSkyrmionThetaDepolTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("theta", "The theta variable");
  params.addParam<Real>("edep", 1.0, "the depolarization field strength");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

EulerSkyrmionThetaDepolTerm::EulerSkyrmionThetaDepolTerm(const InputParameters & parameters)
  :Kernel(parameters),
   _theta_var(coupled("theta")),
   _theta(coupledValue("theta")),
   _edep(getParam<Real>("edep")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
EulerSkyrmionThetaDepolTerm::computeQpResidual()
{
  return _test[_i][_qp] * _edep * std::sin(_theta[_qp]);
}

Real
EulerSkyrmionThetaDepolTerm::computeQpJacobian()
{
  return _test[_i][_qp] * _edep * _phi[_j][_qp] * std::cos(_theta[_qp]);
}

