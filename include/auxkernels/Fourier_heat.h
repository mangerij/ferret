#ifndef Fourier_heat_H
#define Fourier_heat_H

#include "AuxKernel.h"
#include "Material.h"

class Fourier_heat;

template<>
InputParameters validParams<Fourier_heat>();

class Fourier_heat: public AuxKernel
{
public:

Fourier_heat(const InputParameters & parameters);

virtual ~Fourier_heat() {}

protected:
virtual Real computeValue();

private:
  const MaterialProperty<Real> & _thC;
  const unsigned int _component;
  const unsigned int _T_var;
  const VariableValue & _T;
  const VariableGradient & _T_grad;
  const Real _len_scale;

};
#endif
