#ifndef HeatFluxTensor_H
#define HeatFluxTensor_H

#include "AuxKernel.h"
#include "Material.h"
#include "RankTwoTensor.h"

class HeatFluxTensor : public AuxKernel
{
public:
HeatFluxTensor(const InputParameters & parameters);

  static InputParameters validParams();

virtual ~HeatFluxTensor() {}

protected:
virtual Real computeValue();

private:

  const VariableValue & _T;
  const VariableGradient & _T_grad;
  const VariableValue & _potential_E_int;
  const VariableGradient & _potential_E_int_grad;
  const MaterialProperty<RankTwoTensor> & _thC_tensor;//for tensor inclusion
  const MaterialProperty<RankTwoTensor> & _ecC_tensor;
  const MaterialProperty<RankTwoTensor> & _sbC_tensor;
  const unsigned int _component;
  const Real _len_scale;
};

#endif
