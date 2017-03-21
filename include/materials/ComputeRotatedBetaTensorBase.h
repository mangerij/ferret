/**
 * @file   ComputeRotatedBetaTensorBase.h
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

#ifndef COMPUTEROTATEDBETATENSORBASE_H
#define COMPUTEROTATEDBETATENSORBASE_H

#include "ComputeBetaTensorBase.h"
#include "RankTwoTensor.h"

/**
 * ComputeRotatedTensorBase the base class for computing photostrictive tensors
 */
class ComputeRotatedBetaTensorBase : public ComputeBetaTensorBase
{
public:
  ComputeRotatedBetaTensorBase(const InputParameters & parameters);

protected:
  RealVectorValue _Euler_angles;
};

#endif //COMPUTEROTATEDBETATENSORBASE_H
