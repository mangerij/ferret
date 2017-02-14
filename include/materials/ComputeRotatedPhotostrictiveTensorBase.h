
#ifndef COMPUTEROTATEDPHOTOSTRICTIVETENSORBASE_H
#define COMPUTEROTATEDPHOTOSTRICTIVETENSORBASE_H

#include "ComputePhotostrictiveTensorBase.h"

/**
 * ComputeRotatedPhotostrictiveTensorBase is an intermediate base class that rotates the photostrictive tensor based on euler angles.
 */
class ComputeRotatedPhotostrictiveTensorBase : public ComputePhotostrictiveTensorBase
{
public:
  ComputeRotatedPhotostrictiveTensorBase(const InputParameters & parameters);

protected:
  RealVectorValue _Euler_angles;
};

#endif //COMPUTEROTATEDPHOTOSTRICTIVETENSORBASE_H
