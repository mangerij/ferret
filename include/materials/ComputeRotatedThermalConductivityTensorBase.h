#ifndef COMPUTEROTATEDTHERMALCONDUCTIVITYTENSORBASE_H
#define COMPUTEROTATEDTHERMALCONDUCTIVITYTENSORBASE_H

#include "ComputeThermalConductivityTensorBase.h"

class ComputeRotatedThermalConductivityTensorBase;

template <>
InputParameters validParams<ComputeRotatedThermalConductivityTensorBase>();

/**
 * ComputeRotatedThermalConductivityTensorBase is an intermediate base class that rotates the linear
 * thermalconductivity tensor, k_{ij},  based on euler angles.
 */
class ComputeRotatedThermalConductivityTensorBase : public ComputeThermalConductivityTensorBase
{
public:
  ComputeRotatedThermalConductivityTensorBase(const InputParameters & parameters);

protected:
  RealVectorValue _Euler_angles;
};

#endif // COMPUTEROTATEDThermalConductivityTENSORBASE_H
