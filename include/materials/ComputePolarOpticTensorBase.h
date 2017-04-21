

#ifndef COMPUTEPOLAROPTICTENSORBASE_H
#define COMPUTEPOLAROPTICTENSORBASE_H

#include "Material.h"
#include "RankTwoTensor.h"

/**
 * ComputePolarOpticTensorBase the base class for computing polar-optic tensors
 */
class ComputePolarOpticTensorBase : public Material
{
public:
  ComputePolarOpticTensorBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpPolarOpticTensor() = 0;

  std::string _base_name;
  std::string _delta_PO_tensor_name;

  MaterialProperty<RankTwoTensor> & _delta_PO_tensor;

};

#endif //COMPUTEPOLAROPTICTENSORBASE_H
