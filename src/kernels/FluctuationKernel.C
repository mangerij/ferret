/**
 * @file   FluctuationKernel.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Jun 14 12:00:20 2016
 *
 * @brief This Kernel is used to introduce noise in between electric field changes in
 *        a quasi-static hysteresis calculation (see arxiv.org/pdf/1701.02613.pdf)
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
