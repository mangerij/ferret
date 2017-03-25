/**
 * @file   Birefringence.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate the birefringence
 * \delta n = n_o - n_e (unrotated, unstressed)
 *
 */

#include "Birefringence.h"

template<>

InputParameters validParams<Birefringence>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("per1", "first perpendicular direction to propagation");
  params.addRequiredCoupledVar("per2", "second perpendicular direction to propagation");
  return params;
}


Birefringence::Birefringence(const InputParameters & parameters) :
  AuxKernel(parameters),
  _var1(coupledValue("per1")),
  _var2(coupledValue("per2"))
{
}

Real
Birefringence::computeValue()
{
  return _var2[_qp] - _var1[_qp];
}


