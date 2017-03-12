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

template<>

InputParameters validParams<RefractiveIndex>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("index_one", "A rank two tensor is being computed, need two indices");
  params.addRequiredParam<unsigned int>("index_two", "A rank two tensor is being computed, need two indices");
  return params;
}


RefractiveIndex::RefractiveIndex(const InputParameters & parameters) :
  AuxKernel(parameters),
   _index_one(getParam<unsigned int>("index_one")),
   _index_two(getParam<unsigned int>("index_two")),
   _indicatrix_vector(getMaterialProperty<RealVectorValue>("indicatrix")),
   _beta_tensor_ij(getMaterialProperty<RankTwoTensor>("beta_tensor"))
{
}

Real
RefractiveIndex::computeValue()
{
  return _indicatrix_vector[_qp](_index_one) + std::pow(-_indicatrix_vector[_qp](_index_one)*_indicatrix_vector[_qp](_index_two) * _beta_tensor_ij[_qp](_index_one, _index_two),0.5);
}


