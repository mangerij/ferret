#ifndef Electric_flux_H
#define Electric_flux_H

#include "AuxKernel.h"
#include "Material.h"
#include "RankTwoTensor.h"

class Electric_flux;

template <>
InputParameters validParams<Electric_flux>();

class Electric_flux : public AuxKernel
{
public:
  Electric_flux(const InputParameters & parameters);

  virtual ~Electric_flux() {}

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
