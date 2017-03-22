/**
 * @file   Birefringence.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate the birefringence
 * \delta n = n_o - n_e
 *
 */

#include "Birefringence.h"

template<>

InputParameters validParams<Birefringence>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("n_o", "The ordinary axis of the indicatrix");
  params.addRequiredCoupledVar("n_e", "The extraordinary axis of the indicatrix");
  return params;
}


Birefringence::Birefringence(const InputParameters & parameters) :
  AuxKernel(parameters),
  _var1(coupledValue("n_o")),
  _var2(coupledValue("n_e"))
{
}

Real
Birefringence::computeValue()
{
  return _var2[_qp] - _var1[_qp];
}


