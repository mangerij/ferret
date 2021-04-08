#ifndef SeebeckEffect_H
#define SeebeckEffect_H

#include "Kernel.h"

class SeebeckEffect;

template <>
InputParameters validParams<SeebeckEffect>();

class SeebeckEffect : public Kernel
{
public:
  SeebeckEffect(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  // const MaterialProperty<Real> & _ecC;
  const unsigned int _component;
  const unsigned int _potential_E_int_var;
  const VariableValue & _potential_E_int;
  const VariableGradient & _potential_E_int_grad;
  const unsigned int _T_var;
  const VariableValue & _T;
  const VariableGradient & _T_grad;
  const MaterialProperty<Real> & _sbC;
  const Real _len_scale;
};
#endif
