#ifndef BOUNDCHARGE_H
#define BOUNDCHARGE_H

#include "AuxKernel.h"


//Forward declarations
class BoundCharge;

template<>
InputParameters validParams<BoundCharge>();

class BoundCharge : public AuxKernel
{
public:
  BoundCharge(const InputParameters & parameters);

  virtual ~BoundCharge() {}

protected:
  virtual Real computeValue();

private:
  const VariableGradient&  _potential_int_grad;
  const VariableGradient&  _potential_ext_grad;
};

#endif
