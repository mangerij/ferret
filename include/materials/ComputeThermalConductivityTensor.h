#ifndef COMPUTETHERMALCONDUCTIVITYTENSOR_H
#define COMPUTETHERMALCONDUCTIVITYTENSOR_H

#include "RankTwoTensor.h"
#include "ComputeRotatedThermalConductivityTensorBase.h"
#include "libmesh/quadrature.h"

class ComputeThermalConductivityTensor;

template <>
InputParameters validParams<ComputeThermalConductivityTensor>();

/**
 * ComputeThermalConductivityTensor defines a linear ThermalConductivity tensor material object with
 * a given base name.
 */
class ComputeThermalConductivityTensor : public ComputeRotatedThermalConductivityTensorBase
{
public:
  ComputeThermalConductivityTensor(const InputParameters & parameters);

protected:
  virtual void computeQpThermalConductivityTensor();

  /// Individual material information
  RankTwoTensor _kij;
};

#endif // COMPUTETHERMALCONDUCTIVITYTENSOR_H
