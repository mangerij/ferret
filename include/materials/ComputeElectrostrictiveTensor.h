
#ifndef COMPUTEELECTROSTRICTIVETENSOR_H
#define COMPUTEELECTROSTRICTIVETENSOR_H

#include "ComputeElasticityTensor.h"
#include "ComputeRotatedElectrostrictiveTensorBase.h"
#include "ElectrostrictiveTensorR4.h"
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
  ElasticityTensorR4 _Cijkl; //the million dollar question is does this get computed _BEFORE_ ?
  ElectrostrictiveTensorR4 _Qmnkl;
};

#endif //COMPUTEELECTROSTRICTIVETENSOR_H
