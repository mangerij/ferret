#ifndef CHERNSIMONSDENSITY_H
#define CHERNSIMONSDENSITY_H

#include "AuxKernel.h"

//Forward declarations
class ChernSimonsDensity;

template<>
InputParameters validParams<ChernSimonsDensity>();

class ChernSimonsDensity : public AuxKernel
{
public:
  ChernSimonsDensity(const InputParameters & parameters);

  virtual ~ChernSimonsDensity() {}

protected:
  virtual Real computeValue();

private:
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const VariableGradient & _polar_x_grad;
  const VariableGradient & _polar_y_grad;
  const VariableGradient & _polar_z_grad;
};

#endif
