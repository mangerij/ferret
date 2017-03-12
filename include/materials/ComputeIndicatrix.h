/**
 * @file   ComputeIndicatrix.h
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

#ifndef COMPUTEINDICATRIX_H
#define COMPUTEINDICATRIX_H

#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include "ComputeIndicatrixBase.h"
#include "libmesh/quadrature.h"

/**
 * ComputeIndicatrix defines the ellipsoidal indicatrix.
 */
class ComputeIndicatrix : public ComputeIndicatrixBase
{
public:
  ComputeIndicatrix(const InputParameters & parameters);

protected:
  virtual void computeQpIndicatrix();
  const Real _no, _ne;
  RealVectorValue _Euler_angles;
};

#endif //COMPUTEINDICATRIX_H
