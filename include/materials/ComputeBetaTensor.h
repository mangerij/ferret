/**
 * @file   ComputeBetaTensor.h
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


#ifndef COMPUTEBETATENSOR_H
#define COMPUTEBETATENSOR_H

#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include "ComputeBetaTensorBase.h"
#include "libmesh/quadrature.h"

/**
 * ComputeBetaTensor defines an photostrictive tensor material object with a given base name.
 */
class ComputeBetaTensor : public ComputeBetaTensorBase
{
public:
  ComputeBetaTensor(const InputParameters & parameters);

protected:
  virtual void computeQpBetaTensor();
  const MaterialProperty<RankTwoTensor> & _strain;
  const MaterialProperty<RankFourTensor> & _photostrictive_tensor;
};

#endif //COMPUTEPHOTOSTRICTIVETENSOR_H
