/**
 * @file   RefractiveIndex.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate an approximate photoelastic change to the refractive index
 *
 * TODO: adjust
 *
 * where \Delta B_{ij} = p_{ijkl} \varepsilon_{kl}
 *
 * for more information, see Chang (Chp. 12 Handbook of Optics).
 *
 */

#include "RefractiveIndex.h"

template<>

InputParameters validParams<RefractiveIndex>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("component", "component");
  params.addRequiredParam<Real>("n_a", "alpha refractive index");
  params.addRequiredParam<Real>("n_b", "beta refractive index");
  params.addRequiredParam<Real>("n_g", "gamma refractive index");
  params.addRequiredCoupledVar("var1", "the change in this refractive index");
  return params;
}


RefractiveIndex::RefractiveIndex(const InputParameters & parameters) :
  AuxKernel(parameters),
   _component(getParam<unsigned int>("component")),
   _na(getParam<Real>("n_a")),
   _nb(getParam<Real>("n_b")),
   _ng(getParam<Real>("n_g")),
  _var1(coupledValue("var1"))
{
}

Real
RefractiveIndex::computeValue()
{
  // the diagonals are related to the B1, B2, B3 terms in rotated indicatrix
//std::pow(  (1.0 / ( _beta_tensor[_qp](_index_i, _index_j)  ) ), 3.0)
  if (_component == 0)
    return _na - _var1[_qp];
  else if (_component == 1)
    return _nb - _var1[_qp];
  else if (_component == 2)
    return _ng - _var1[_qp];
  else
    return 0.0;
}


