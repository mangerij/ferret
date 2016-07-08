/**
 * @file   KarmanenkoDriver.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * named after Karmanenko et al J. Euro. Ceram. Soc. 27 (2007) 3109–3112
 * this term drives the temperature changes due to the field-induced entropic changes
 * NOTE: this is just a test kernel for now, as the anisotropy of _grad_potential_int
 * needs to be taken into account
 *
 * The procedure is as follows, dEstep will be related to the stepping procedure in the
 * quasi-static hysteresis loop. The only difficulty will pinning down how noise 
 * introduced is related to this kernel.
 *
 * Currently the kernel is setup for just a field along z and using the approximation of 
 * Gu et al Appl. Phys. Lett. 102, 112901, (2013).
*/


#ifndef KARMANENKODRIVER_H
#define KARMANENKODRIVER_H

#include "Kernel.h"

class KarmanenkoDriver;

template<>
InputParameters validParams<KarmanenkoDriver>();

class KarmanenkoDriver: public Kernel
{
public:

  KarmanenkoDriver(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _potential_int_var;
  const unsigned int _potential_ext_var;
  const VariableGradient &  _potential_int_grad;
  const VariableGradient &  _potential_ext_grad;
  const unsigned int _temperature_var;
  const VariableValue & _temperature;
  const Real _rho1, _C1, _C2, _C3, _C4, _dEstep;
  const Real _len_scale;
};
#endif //KARMANENKODRIVER_H
