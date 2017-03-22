#ifndef POLARMAG_H
#define POLARMAG_H

#include "AuxKernel.h"

//Forward Declarations
class PolarMag;

template<>
InputParameters validParams<PolarMag>();

/**
 * Coupled auxiliary value
 */
class PolarMag: public AuxKernel
{
public:
  PolarMag(const InputParameters & parameters);

  virtual ~PolarMag() {}

protected:
    virtual Real computeValue();

private:
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
};

#endif // NEIGHBORANGLEMARKER_H
