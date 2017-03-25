/**
 * @file   RefractiveIndex.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate an approximate photoelastic change to the refractive index
 * \Delta \epsilon_{ij} = - n_i^2 n_j^2 \Delta B_{ij}
 *
 * where \Delta B_{ij} = p_{ijkl} \varepsilon_{kl}
 *
 * for more information, see Chang (Chp. 12 Handbook of Optics).
 *
 */

#include "RefractiveIndex.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

template<>

InputParameters validParams<RefractiveIndex>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("index_i", "index of the delta beta vector");
  params.addRequiredParam<Real>("n", "need index of refraction along this (i, j) principle direction");
  return params;
}


RefractiveIndex::RefractiveIndex(const InputParameters & parameters) :
  AuxKernel(parameters),
   _index_i(getParam<unsigned int>("index_i")),
   _n(getParam<Real>("n")),
   _delta_beta_tensor(getMaterialProperty<RealTensorValue>("delta_beta_tensor"))
{
}

Real
RefractiveIndex::computeValue()
{
  return - 0.5 * std::pow(_n, 3.0) * _delta_beta_tensor[_qp](0, _index_i); //note n needs to be the unstressed indicatrix entry. It is still unclear what B4, B5, B6 provide to this.
}


