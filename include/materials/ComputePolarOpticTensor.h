/**
 * @file   ComputePolarOpticTensor.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * Calculate an approximate polar-optic change to the refractive index
 * \delta B_{ij} = p_{ijkl} Q_{klmn} P_m P_n
 *
 */

#ifndef COMPUTEPOLAROPTICTENSOR_H
#define COMPUTEPOLAROPTICTENSOR_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "ComputePolarOpticTensorBase.h"

/**
 * ComputePolarOpticTensor the base class for computing polar-optic adjustments to B_{ij}
 */
class ComputePolarOpticTensor : public ComputePolarOpticTensorBase
{
public:
  ComputePolarOpticTensor(const InputParameters & parameters);

protected:
  virtual void computeQpPolarOpticTensor();

  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;

  const MaterialProperty<RankFourTensor> & _photostrictive_tensor;
  const MaterialProperty<RankFourTensor> & _electrostrictive_tensor;
};

#endif //COMPUTEPOLAROPTICTENSOR_H
