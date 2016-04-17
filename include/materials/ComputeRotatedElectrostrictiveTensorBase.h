
#ifndef COMPUTEROTATEDELECTROSTRICTIVETENSORBASE_H
#define COMPUTEROTATEDELECTROSTRICTIVETENSORBASE_H

#include "ComputeElectrostrictiveTensorBase.h"

/**
 * ComputeRotatedElectrostrictiveTensorBase is an intermediate base class that rotates the electrostrictive tensor based on euler angles.
 */
class ComputeRotatedElectrostrictiveTensorBase : public ComputeElectrostrictiveTensorBase
{
public:
  ComputeRotatedElectrostrictiveTensorBase(const InputParameters & parameters);

protected:
  RealVectorValue _Euler_angles;
};

#endif //COMPUTEROTATEDELECTROSTRICTIVETENSORBASE_H
