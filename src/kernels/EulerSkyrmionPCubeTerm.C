/**
 * @file   EulerSkyrmionPCubeTerm.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Feb. 20. 2017
 *
 */

#include "EulerSkyrmionPCubeTerm.h"

class EulerSkyrmionPCubeTerm;

template<>
InputParameters validParams<EulerSkyrmionPCubeTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("P", "The polar magnitude variable");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("P0", 1.0, "the domain wall coefficient");
  return params;
}

EulerSkyrmionPCubeTerm::EulerSkyrmionPCubeTerm(const InputParameters & parameters)
  :Kernel(parameters),
   _P_var(coupled("P")),
   _P(coupledValue("P")),
   _len_scale(getParam<Real>("len_scale")),
   _P0(getParam<Real>("P0"))
{
}

Real
EulerSkyrmionPCubeTerm::computeQpResidual()
{
  return _test[_i][_qp] * _P[_qp] * _P[_qp] * _P[_qp] / (_P0 * _P0);
}

Real
EulerSkyrmionPCubeTerm::computeQpJacobian()
{
  return 2.0 * _test[_i][_qp] * _phi[_j][_qp] * _P[_qp] * _P[_qp] / (_P0 * _P0);
}
