#ifndef RESIDUALV_TENSOR_H
#define RESIDUALV_TENSOR_H

#include "Kernel.h"
#include "Material.h"
#include "RankTwoTensor.h" //added for tensor calculation

class ResidualV_tensor;

template <>
InputParameters validParams<ResidualV_tensor>();

class ResidualV_tensor : public Kernel
{
public:
  ResidualV_tensor(const InputParameters & parameters);

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
  const MaterialProperty<RankTwoTensor> & _ecC_tensor; // for tensor inclusion
  const MaterialProperty<RankTwoTensor> & _sbC_tensor;
  const Real _len_scale;
};
#endif
