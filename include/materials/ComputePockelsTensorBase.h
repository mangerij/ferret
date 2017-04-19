#ifndef COMPUTEPOCKELSTENSORBASE_H
#define COMPUTEPOCKELSTENSORBASE_H

#include "Material.h"
#include "RankFourTensor.h"

/**
 * ComputePockelsTensorBase the base class for computing photostrictive tensors
 */
class ComputePockelsTensorBase : public Material
{
public:
  ComputePockelsTensorBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpPockelsTensor() = 0;

  std::string _base_name;
  std::string _pockels_tensor_name;

  MaterialProperty<RankFourTensor> & _pockels_tensor;
};

#endif //COMPUTEPHOTOSTRICTIVETENSORBASE_H
