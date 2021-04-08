#ifndef zT_aux_H
#define zT_aux_H

#include "AuxKernel.h"


//Forward declarations
class zT_aux;

template<>
InputParameters validParams<zT_aux>();


class zT_aux : public AuxKernel
{
public:
  zT_aux(const InputParameters & parameters);

  virtual ~zT_aux() {}

protected:
  virtual Real computeValue();

private:
  const unsigned int _T_var;
  const VariableValue & _T;
  const VariableGradient & _T_grad;
  // const Real _electrical_conductivity;
  // const Real _seebeck_coefficient;
  // const Real _thermal_conductivity;
  const MaterialProperty<Real> & _ecC;
  const MaterialProperty<Real> & _sbC;
  const MaterialProperty<Real> & _thC;
};
#endif
