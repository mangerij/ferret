/**
 * @file   PrefactorRefractive.C
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

#include "PrefactorRefractive.h"
#include "RotationTensor.h"
#include "RankTwoTensor.h"

template<>

InputParameters validParams<PrefactorRefractive>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("index_i", "first index of the beta vector");
  params.addRequiredParam<unsigned int>("index_j", "second index of the beta vector");
  params.addRequiredParam<unsigned int>("index_k", "first index of the delta_beta vector");
  params.addRequiredParam<unsigned int>("index_l", "second index of the delta_beta vector");
  return params;
}


PrefactorRefractive::PrefactorRefractive(const InputParameters & parameters) :
  AuxKernel(parameters),
   _index_i(getParam<unsigned int>("index_i")),
   _index_j(getParam<unsigned int>("index_j")),
   _index_k(getParam<unsigned int>("index_k")),
   _index_l(getParam<unsigned int>("index_l")),
   _beta_tensor(getMaterialProperty<RankTwoTensor>("beta_tensor")),
   _delta_beta_tensor(getMaterialProperty<RealTensorValue>("delta_beta_tensor")) //or need RealTensorValue
{
}

Real
PrefactorRefractive::computeValue()
{
  // the diagonals are related to the B1, B2, B3 terms in rotated indicatrix
//std::pow(  (1.0 / ( _beta_tensor[_qp](_index_i, _index_j)  ) ), 3.0) 
  return _beta_tensor[_qp](_index_k, _index_l);
//- 0.5 * std::pow((1.0 / std::pow(std::abs(_beta_tensor[_qp](_index_k, _index_l)), 0.5), 3.0);
}


