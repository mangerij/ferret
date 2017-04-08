/**
 * @file   ComputeBetaTensorBase.h
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

#ifndef COMPUTEBETATENSORBASE_H
#define COMPUTEBETATENSORBASE_H

#include "Material.h"
#include "RankTwoTensor.h"

/**
 * ComputeTensorBase the base class for computing photostrictive tensors
 */
class ComputeBetaTensorBase : public Material
{
public:
  ComputeBetaTensorBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpBetaTensor() = 0;

  std::string _base_name;
  std::string _beta_tensor_name;

  MaterialProperty<RankTwoTensor> & _beta_tensor;

};

#endif //COMPUTEBETATENSORBASE_H
