/**
 * @file   PrefactorRefractive.h
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

#ifndef PREFACTORREFRACTIVE_H
#define PREFACTORREFRACTIVE_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class PrefactorRefractive;

template<>
InputParameters validParams<PrefactorRefractive>();


class PrefactorRefractive : public AuxKernel
{
public:
  PrefactorRefractive(const InputParameters & parameters);

  virtual ~PrefactorRefractive() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _index_i;
  const unsigned int _index_j;
  const MaterialProperty<RankTwoTensor> & _beta_tensor;
};

#endif // CHANGEINREFRACTIVEINDEX_H
