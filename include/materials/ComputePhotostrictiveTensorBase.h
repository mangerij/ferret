
#ifndef COMPUTEPHOTOSTRICTIVETENSORBASE_H
#define COMPUTEPHOTOSTRICTIVETENSORBASE_H

#include "Material.h"
#include "RankFourTensor.h"

/**
 * ComputePhotostrictiveTensorBase the base class for computing photostrictive tensors
 */
class ComputePhotostrictiveTensorBase : public Material
{
public:
  ComputePhotostrictiveTensorBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpPhotostrictiveTensor() = 0;

  std::string _base_name;
  std::string _photostrictive_tensor_name;
  std::string _photostrictive_tensorP_name;

  MaterialProperty<RankFourTensor> & _photostrictive_tensor;
  MaterialProperty<RankFourTensor> & _photostrictive_tensorP;
};

#endif //COMPUTEPHOTOSTRICTIVETENSORBASE_H
