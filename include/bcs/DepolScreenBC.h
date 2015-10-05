

#ifndef DEPOLSCREENBC_H
#define DEPOLSCREENBC_H

#include "NodalNormalBC.h"

class DepolScreenBC;

template<>
InputParameters validParams<DepolScreenBC>();

/**
 * Boundary condition of a Dirichlet type but ON the P dot n value.
 * Sets the value in the node based on the screening from the electrode on the
 * unscreened depolarization field
 */

class DepolScreenBC : public NodalNormalBC
{
public:
  DepolScreenBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const Real & _value;

};

#endif /* DEPOLSCREENBC_H */
