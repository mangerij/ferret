/**
 * @file   ComputeDeltaBetaTensor.h
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


#ifndef COMPUTEDELTABETATENSOR_H
#define COMPUTEDELTABETATENSOR_H

#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include "ComputeDeltaBetaTensorBase.h"
#include "libmesh/quadrature.h"

/**
 * ComputeDeltaBetaTensor defines an impermeability tensor material object with a given base name.
 */
class ComputeDeltaBetaTensor : public ComputeDeltaBetaTensorBase
{
public:
  ComputeDeltaBetaTensor(const InputParameters & parameters);

protected:
  virtual void computeQpDeltaBetaTensor();
  const MaterialProperty<RankTwoTensor> & _strain;
  const MaterialProperty<RankFourTensor> & _photostrictive_tensor;
};

#endif //COMPUTEDELTABETATENSOR_H
