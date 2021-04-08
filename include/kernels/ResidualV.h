#ifndef RESIDUALV_H
#define RESIDUALV_H

#include "Kernel.h"

class ResidualV;

template <>
InputParameters validParams<ResidualV>();

class ResidualV : public Kernel
{
public:
  ResidualV(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _potential_E_int_var;
  const VariableValue & _potential_E_int;
  const VariableGradient & _potential_E_int_grad;
  const unsigned int _T_var;
  const VariableValue & _T;
  const VariableGradient & _T_grad;
  const MaterialProperty<Real> & _ecC;
  const MaterialProperty<Real> & _sbC;
};
#endif
