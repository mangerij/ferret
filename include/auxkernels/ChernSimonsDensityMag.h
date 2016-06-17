#ifndef CHERNSIMONSDENSITYMAG_H
#define CHERNSIMONSDENSITYMAG_H

#include "AuxKernel.h"

//Forward declarations
class ChernSimonsDensityMag;

template<>
InputParameters validParams<ChernSimonsDensityMag>();

class ChernSimonsDensityMag : public AuxKernel
{
public:
  ChernSimonsDensityMag(const InputParameters & parameters);

  virtual ~ChernSimonsDensityMag() {}

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
