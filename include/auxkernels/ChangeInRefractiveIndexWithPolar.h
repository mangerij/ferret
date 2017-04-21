/**
 * @file   ChangeInRefractiveIndexWithPolar.h
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

#ifndef CHANGEINREFRACTIVEINDEXWITHPOLAR_H
#define CHANGEINREFRACTIVEINDEXWITHPOLAR_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class ChangeInRefractiveIndexWithPolar;

template<>
InputParameters validParams<ChangeInRefractiveIndexWithPolar>();


class ChangeInRefractiveIndexWithPolar : public AuxKernel
{
public:
  ChangeInRefractiveIndexWithPolar(const InputParameters & parameters);

  virtual ~ChangeInRefractiveIndexWithPolar() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _index_i;
  const unsigned int _index_j;
  const unsigned int _index_k;
  const unsigned int _index_l;
  const MaterialProperty<RankTwoTensor> & _beta_tensor;
  const MaterialProperty<RankTwoTensor> & _delta_beta_tensor;
  const MaterialProperty<RankTwoTensor> & _delta_PO_tensor;
};

#endif // CHANGEINREFRACTIVEINDEXWITHPOLAR_H
