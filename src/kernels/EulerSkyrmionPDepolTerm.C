/**
 * @file   EulerSkyrmionPDepolTerm.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Feb. 20. 2017
 *
 */

#include "EulerSkyrmionPDepolTerm.h"

class EulerSkyrmionPDepolTerm;

template<>
InputParameters validParams<EulerSkyrmionPDepolTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("theta", "The theta variable");
  params.addParam<Real>("edep", 1.0, "the depolarization field strength");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

EulerSkyrmionPDepolTerm::EulerSkyrmionPDepolTerm(const InputParameters & parameters)
  :Kernel(parameters),
   _theta_var(coupled("theta")),
   _theta(coupledValue("theta")),
   _edep(getParam<Real>("edep")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
EulerSkyrmionPDepolTerm::computeQpResidual()
{
  return - _test[_i][_qp] * _edep * std::cos(_theta[_qp]);
}

Real
EulerSkyrmionPDepolTerm::computeQpJacobian()
{
  return 0.0;
}

Real
EulerSkyrmionPDepolTerm::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _theta_var)
    return _test[_i][_qp] * _edep * std::sin(_theta[_qp]);
  else
    return 0.0;
}
