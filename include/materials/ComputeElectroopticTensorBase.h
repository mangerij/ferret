#ifndef COMPUTEELECTROOPTICTENSORBASE_H
#define COMPUTEELECTROOPTICTENSORBASE_H

#include "Material.h"
#include "RankThreeTensor.h"

/**
 * ComputeElectroopticTensorBase the base class for computing photostrictive tensors
 */
class ComputeElectroopticTensorBase : public Material
{
public:
  ComputeElectroopticTensorBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpElectroopticTensor() = 0;

  std::string _base_name;
  std::string _electrooptic_tensor_name;

  MaterialProperty<RankThreeTensor> & _electrooptic_tensor;
};

#endif //COMPUTEELECTROOPTICTENSORBASE_H
