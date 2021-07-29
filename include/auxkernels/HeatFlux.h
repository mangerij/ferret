#ifndef HeatFlux_H
#define HeatFlux_H

#include "AuxKernel.h"
#include "Material.h"
#include "RankTwoTensor.h"

class HeatFlux;

template <>
InputParameters validParams<HeatFlux>();

class HeatFlux : public AuxKernel
{
public:
  HeatFlux(const InputParameters & parameters);

  virtual ~HeatFlux() {}

protected:
  virtual Real computeValue();

private:
  const VariableValue & _T;
  const VariableGradient & _T_grad;
  const VariableValue & _potential_E_int;
  const VariableGradient & _potential_E_int_grad;
  const MaterialProperty<Real> & _thC;
  const MaterialProperty<Real> & _ecC;
  const MaterialProperty<Real> & _sbC;
  const unsigned int _component;
};

#endif
