/**
 * @file   RefractiveIndex.h
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

#ifndef REFRACTIVEINDEX_H
#define REFRACTIVEINDEX_H

#include "AuxKernel.h"

//Forward declarations
class RefractiveIndex;

template<>
InputParameters validParams<RefractiveIndex>();


class RefractiveIndex : public AuxKernel
{
public:
  RefractiveIndex(const InputParameters & parameters);

  virtual ~RefractiveIndex() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _component;
  Real _na;
  Real _nb;
  Real _ng;
  const VariableValue & _var1;
};

#endif // CHANGEINREFRACTIVEINDEX_H
