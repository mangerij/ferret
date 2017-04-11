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
  return params;
}


PrefactorRefractive::PrefactorRefractive(const InputParameters & parameters) :
  AuxKernel(parameters),
   _index_i(getParam<unsigned int>("index_i")),
   _index_j(getParam<unsigned int>("index_j")),
   _beta_tensor(getMaterialProperty<RankTwoTensor>("beta_tensor"))
{
}

Real
PrefactorRefractive::computeValue()
{
  //Moose::out << "\n B"; std::cout << _index_i; std::cout << _index_j; Moose::out << " = "; std::cout << _beta_tensor[_qp](_index_i,_index_j);
  return _beta_tensor[_qp](_index_i, _index_j);
}


