/**
 * @file   ComputeIndicatrixBase.h
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


#ifndef COMPUTEINDICATRIXBASE_H
#define COMPUTEINDICATRIXBASE_H

#include "Material.h"
#include "RankTwoTensor.h"

/**
 * ComputeIndicatrixBase the base class for computing the indicatrix.
 */
class ComputeIndicatrixBase : public Material
{
public:
  ComputeIndicatrixBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpIndicatrix() = 0;

  std::string _base_name;
  std::string _indicatrix_name;

  MaterialProperty<RealVectorValue> & _indicatrix;
};

#endif //COMPUTEINDICATRIXBASE_H
