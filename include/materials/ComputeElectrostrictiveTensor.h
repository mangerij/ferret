
#ifndef COMPUTEELECTROSTRICTIVETENSOR_H
#define COMPUTEELECTROSTRICTIVETENSOR_H

#include "RankFourTensor.h"
// #include "ComputeElasticityTensor.h"
#include "ElectrostrictiveTensorTools.h"
#include "ComputeRotatedElectrostrictiveTensorBase.h"
#include "libmesh/quadrature.h"

/**
 * ComputeElectrostrictiveTensor defines an electrostrictive tensor material object with a given base name.
 */
class ComputeElectrostrictiveTensor : public ComputeRotatedElectrostrictiveTensorBase
{
public:
  ComputeElectrostrictiveTensor(const InputParameters & parameters);

protected:
  virtual void computeQpElectrostrictiveTensor();

  /// Individual material information
  RankFourTensor _Qmnkl;
  RankFourTensor _qijkl;

 private:
   const MaterialProperty<RankFourTensor> & _elasticity_tensor;
};

#endif //COMPUTEELECTROSTRICTIVETENSOR_H
