#ifndef COMPUTEELECTROOPTICTENSOR_H
#define COMPUTEELECTROOPTICTENSOR_H

#include "RankThreeTensor.h"
#include "ComputeRotatedElectroopticTensorBase.h"

/**
 * ComputeElectroopticTensor defines a linear electrooptic tensor material object with a given base name.
 */
class ComputeElectroopticTensor : public ComputeRotatedElectroopticTensorBase
{
public:
  ComputeElectroopticTensor(const InputParameters & parameters);

protected:
  virtual void computeQpElectroopticTensor();

  /// Individual material information
  RankThreeTensor _rijk;
};

#endif //COMPUTEELECTROOPTICTENSOR_H
