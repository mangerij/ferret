#ifndef COMPUTEPHOTOSTRICTIVETENSOR_H
#define COMPUTEPHOTOSTRICTIVETENSOR_H

#include "RankFourTensor.h"
#include "ComputeRotatedPhotostrictiveTensorBase.h"
#include "libmesh/quadrature.h"

/**
 * ComputePhotostrictiveTensor defines an photostrictive tensor material object with a given base name.
 */
class ComputePhotostrictiveTensor : public ComputeRotatedPhotostrictiveTensorBase
{
public:
  ComputePhotostrictiveTensor(const InputParameters & parameters);

protected:
  virtual void computeQpPhotostrictiveTensor();

  /// Individual material information
  RankFourTensor _Pmnkl;
};

#endif //COMPUTEPHOTOSTRICTIVETENSOR_H
