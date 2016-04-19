
#ifndef COMPUTEELECTROSTRICTIVETENSORBASE_H
#define COMPUTEELECTROSTRICTIVETENSORBASE_H

#include "Material.h"
#include "RankFourTensor.h"

/**
 * ComputeElectrostrictiveTensorBase the base class for computing electrostrictive tensors
 */
class ComputeElectrostrictiveTensorBase : public Material
{
public:
  ComputeElectrostrictiveTensorBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpElectrostrictiveTensor() = 0;

  std::string _base_name;
  std::string _electrostrictive_tensor_name;

  MaterialProperty<RankFourTensor> & _electrostrictive_tensor;
};

#endif //COMPUTEELECTROSTRICTIVETENSORBASE_H
