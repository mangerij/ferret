/**
 * @file   Impermeability.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * for more information, see Chang (Chp. 12 Handbook of Optics).
 *
 */

#include "Impermeability.h"
#include "RankTwoTensor.h"

template<>

InputParameters validParams<Impermeability>()

{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("index_i", "A rank two tensor is being computed, need two indices");
  params.addRequiredParam<unsigned int>("index_j", "A rank two tensor is being computed, need two indices");
  return params;
}


Impermeability::Impermeability(const InputParameters & parameters) :
  AuxKernel(parameters),
   _index_i(getParam<unsigned int>("index_i")),
   _index_j(getParam<unsigned int>("index_j")),
   _delta_beta_tensor(getMaterialProperty<RankTwoTensor>("delta_beta_tensor")),
   _beta_tensor(getMaterialProperty<RankTwoTensor>("beta_tensor"))
{
}

Real
Impermeability::computeValue()
{
  return _beta_tensor[_qp](_index_i, _index_j) + _delta_beta_tensor[_qp](_index_i, _index_j);
}


