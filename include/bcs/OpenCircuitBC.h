/****************************************************************/
/* Hydrostatics BC:                                             */
/*     This BC is intended to implement                         */
/*     P + \epsilon * _grad\phi = 0 at a boundary               */
/****************************************************************/

#ifndef OPENCIRCUITBC_H
#define OPENCIRCUITBC_H

#include "IntegratedBC.h"

//Forward Declarations
class OpenCircuitBC;

template<>
InputParameters validParams<OpenCircuitBC>();

class OpenCircuitBC : public IntegratedBC
{
public:

  OpenCircuitBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

private:
  const unsigned int _component;
  const VariableGradient &  _potential_int_grad;
  const Real _permittivity;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
};

#endif
