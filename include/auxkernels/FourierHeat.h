#ifndef FourierHeat_H
#define FourierHeat_H

#include "AuxKernel.h"
#include "Material.h"

class FourierHeat;

class FourierHeat: public AuxKernel
{
public:

FourierHeat(const InputParameters & parameters);

static InputParameters validParams();

virtual ~FourierHeat() {}

protected:
virtual Real computeValue();

private:
  const MaterialProperty<Real> & _thC;
  const unsigned int _component;
  const unsigned int _T_var;
  const VariableValue & _T;
  const VariableGradient & _T_grad;
  const Real _len_scale;

};
#endif
