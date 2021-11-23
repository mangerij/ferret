#ifndef ELECTRICFLUX_H
#define ELECTRICFLUX_H

#include "AuxKernel.h"
#include "Material.h"
#include "RankTwoTensor.h"

class ElectricFlux : public AuxKernel
{
public:
  ElectricFlux(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~ElectricFlux() {}

protected:
  virtual Real computeValue();

private:
  const VariableValue & _T;
  const VariableGradient & _T_grad;
  const VariableValue & _potential_E_int;
  const VariableGradient & _potential_E_int_grad;
  const MaterialProperty<Real> & _ecC;
  const MaterialProperty<Real> & _sbC;
  const unsigned int _component;
};

#endif
