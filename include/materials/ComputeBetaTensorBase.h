#ifndef COMPUTEBETATENSORBASE_H
#define COMPUTEBETATENSORBASE_H

#include "Material.h"
#include "RankTwoTensor.h"

/**
 * ComputeBetaTensorBase the base class for computing photostrictive tensors
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
