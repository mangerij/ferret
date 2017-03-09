#ifndef DIVP_H
#define DIVP_H

#include "AuxKernel.h"

//Forward declarations
class DivP;

template<>
InputParameters validParams<DivP>();

class DivP : public AuxKernel
{
public:
  DivP(const InputParameters & parameters);

  virtual ~DivP() {}

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
