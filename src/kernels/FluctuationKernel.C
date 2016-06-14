/**
 * @file   FluctuationKernel.C
 * @author J. Mangeri <john.mangeri@uconn.edu
 *
 */

#include "FluctuationKernel.h"
#include<cmath>

template<>
InputParameters validParams<FluctuationKernel>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("deltaPi", 0.0, "The magnitude of the fluctuation across the ith component");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  return params;
}

FluctuationKernel::FluctuationKernel(const InputParameters & parameters)
  :Kernel(parameters),
   _deltaPi(coupledValue("deltaPi")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
FluctuationKernel::computeQpResidual()
{
  return -_deltaPi[_qp] * _test[_i][_qp];
}

Real
FluctuationKernel::computeQpJacobian()
{
 
  return 0.0;
}

Real
FluctuationKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
    return 0.0;
}
