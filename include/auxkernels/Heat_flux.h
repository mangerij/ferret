#ifndef Heat_flux_H
#define Heat_flux_H

#include "AuxKernel.h"
#include "Material.h"
#include "RankTwoTensor.h"

class Heat_flux;

template <>
InputParameters validParams<Heat_flux>();

class Heat_flux : public AuxKernel
{
public:
  Heat_flux(const InputParameters & parameters);

  virtual ~Heat_flux() {}

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
  const Real _len_scale;
};

#endif
