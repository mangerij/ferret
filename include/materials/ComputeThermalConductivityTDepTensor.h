#ifndef COMPUTETHERMALCONDUCTIVITYTDEPTENSOR_H
#define COMPUTETHERMALCONDUCTIVITYTDEPTENSOR_H

#include "RankTwoTensor.h"
#include "ComputeRotatedThermalConductivityTensorBase.h"
#include "libmesh/quadrature.h"

/**
 * ComputeThermalConductivityTensor defines a linear ThermalConductivity tensor material object with
 * a given base name.
 */

class ComputeThermalConductivityTDepTensor : public ComputeRotatedThermalConductivityTensorBase
{
public:
  ComputeThermalConductivityTDepTensor(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpThermalConductivityTensor();

  /// Individual material information
  const VariableValue & _T;

  RankTwoTensor _akij;
  RankTwoTensor _bkij;
  RankTwoTensor _ckij;
};

#endif // COMPUTETHERMALCONDUCTIVITYTDEPTENSOR_H
