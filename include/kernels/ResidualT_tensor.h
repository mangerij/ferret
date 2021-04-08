#ifndef RESIDUALT_TENSOR_H
#define RESIDUALT_TENSOR_H

#include "Kernel.h"
#include "Material.h"
#include "RankTwoTensor.h" //added for tensor calculation

class ResidualT_tensor;

template <>
InputParameters validParams<ResidualT_tensor>();

class ResidualT_tensor : public Kernel
{
public:
  ResidualT_tensor(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _component;
  const unsigned int _potential_E_int_var;
  const VariableValue & _potential_E_int;
  const VariableGradient & _potential_E_int_grad;
  const unsigned int _T_var;
  const VariableValue & _T;
  const VariableGradient & _T_grad;
  const MaterialProperty<RankTwoTensor> & _thC_tensor; // for tensor inclusion
  const MaterialProperty<RankTwoTensor> & _ecC_tensor;
  const MaterialProperty<RankTwoTensor> & _sbC_tensor;
  const Real _len_scale;
};
#endif
