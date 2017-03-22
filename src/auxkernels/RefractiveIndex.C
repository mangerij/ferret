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
  params.addRequiredParam<unsigned int>("index_one", "A rank two tensor is being computed, need three indices");
  params.addRequiredParam<unsigned int>("index_two", "A rank two tensor is being computed, need three indices");
  params.addRequiredParam<unsigned int>("index_three", "A rank two tensor is being computed, need three indices");
  return params;
}


RefractiveIndex::RefractiveIndex(const InputParameters & parameters) :
  AuxKernel(parameters),
   _index_one(getParam<unsigned int>("index_one")),
   _index_two(getParam<unsigned int>("index_two")),
   _index_three(getParam<unsigned int>("index_three")),
   _indicatrix_vector(getMaterialProperty<RealVectorValue>("indicatrix")),
   _beta_tensor_ij(getMaterialProperty<RankTwoTensor>("beta_tensor"))
{
}

Real
RefractiveIndex::computeValue()
{
  return (_indicatrix_vector[_qp](_index_one) 
+ _indicatrix_vector[_qp](_index_two) 
+ _indicatrix_vector[_qp](_index_three) )/3.0; //+ std::pow(-std::pow(_indicatrix_vector[_qp](_index_one), 2) * std::pow(_indicatrix_vector[_qp](_index_two), 2) * _beta_tensor_ij[_qp](_index_one, _index_one), 0.5);
}


