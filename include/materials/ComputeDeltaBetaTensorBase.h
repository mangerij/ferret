/**
 * @file   ComputeDeltaBetaTensorBase.h
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

#ifndef COMPUTEDELTABETATENSORBASE_H
#define COMPUTEDELTABETATENSORBASE_H

#include "Material.h"
#include "RankTwoTensor.h"

/**
 * ComputeBetaTensorBase the base class for computing photostrictive tensors
 */
class ComputeDeltaBetaTensorBase : public Material
{
public:
  ComputeDeltaBetaTensorBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpDeltaBetaTensor() = 0;

  std::string _base_name;
  std::string _delta_beta_tensor_name;

  MaterialProperty<RankTwoTensor> & _delta_beta_tensor;

};

#endif //COMPUTEDELTABETATENSORBASE_H
