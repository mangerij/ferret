/**
 * @file   ChangeInRefractiveIndex.h
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

#ifndef CHANGEINREFRACTIVEINDEX_H
#define CHANGEINREFRACTIVEINDEX_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class ChangeInRefractiveIndex;

template<>
InputParameters validParams<ChangeInRefractiveIndex>();


class ChangeInRefractiveIndex : public AuxKernel
{
public:
  ChangeInRefractiveIndex(const InputParameters & parameters);

  virtual ~ChangeInRefractiveIndex() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _index_i;
  const unsigned int _index_j;
  const unsigned int _index_k;
  const unsigned int _index_l;
  const MaterialProperty<RankTwoTensor> & _beta_tensor;
  const MaterialProperty<RankTwoTensor> & _delta_beta_tensor;
};

#endif // CHANGEINREFRACTIVEINDEX_H
