
#ifndef COMPUTEROTATEDELECTROOPTICTENSORBASE_H
#define COMPUTEROTATEDELECTROOPTICTENSORBASE_H

#include "ComputeElectroopticTensorBase.h"

/**
 * ComputeRotatedElectroopticTensorBase is an intermediate base class that rotates the linear electrooptic tensor, r_{ijk},  based on euler angles.
 */
class ComputeRotatedElectroopticTensorBase : public ComputeElectroopticTensorBase
{
public:
  ComputeRotatedElectroopticTensorBase(const InputParameters & parameters);

protected:
  RealVectorValue _Euler_angles;
};

#endif //COMPUTEROTATEDELECTROOPTICTENSORBASE_H
