#ifndef PKNORM_H
#define PKNORM_H

#include "AuxKernel.h"

//Forward Declarations
class PkNorm;

template<>
InputParameters validParams<PkNorm>();

/**
 * Coupled auxiliary value
 */
class PkNorm: public AuxKernel
{
public:
  PkNorm(const InputParameters & parameters);

  virtual ~PkNorm() {}

protected:
    virtual Real computeValue();

private:
  const unsigned int _component;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
};

#endif // PKNORM_H
