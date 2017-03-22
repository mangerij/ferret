/**
 * @file   ComputeBetaTensor.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate an approximate photoelastic change to the refractive index
 *
 * where \Delta (1/n^2) = \Delta B_{ij} = p_{ijkl} \varepsilon_{kl}
 *
 * Note that B_{ij} = \epsilon_{ij}^{-1} in the principle axis frame.
 *
 * for more information, see Chang (Chp. 12 Handbook of Optics).
 *
 */

#ifndef COMPUTEBETATENSOR_H
#define COMPUTEBETATENSOR_H

#include "RankTwoTensor.h"
#include "ComputeRotatedBetaTensorBase.h"
#include "libmesh/quadrature.h"

/**
 * ComputeDeltaBetaTensor defines an impermeability tensor material object with a given base name.
 */
class ComputeBetaTensor : public ComputeRotatedBetaTensorBase
{
public:
  ComputeBetaTensor(const InputParameters & parameters);

protected:
  virtual void computeQpBetaTensor();
  /// Individual material information
  Real _no;
  Real _ne;
};

#endif //COMPUTEBETATENSOR_H
