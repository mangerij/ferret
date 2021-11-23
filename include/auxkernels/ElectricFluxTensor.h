#ifndef ElectricFluxTensor_H
#define ElectricFluxTensor_H

#include "AuxKernel.h"
#include "Material.h"
#include "RankTwoTensor.h"

class ElectricFluxTensor : public AuxKernel
{
public:
ElectricFluxTensor(const InputParameters & parameters);

  static InputParameters validParams();

virtual ~ElectricFluxTensor() {}

protected:
virtual Real computeValue();

private:

  const VariableValue & _T;
  const VariableGradient & _T_grad;
  const VariableValue & _potential_E_int;
  const VariableGradient & _potential_E_int_grad;
  const MaterialProperty<RankTwoTensor> & _ecC_tensor;//tensor inclusion
  const MaterialProperty<RankTwoTensor> & _sbC_tensor;
  const unsigned int _component;
};

#endif
