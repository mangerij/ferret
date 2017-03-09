/**
 * @file   AnisotropicElectrostatics.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 */

#ifndef ANISOTROPICELECTROSTATICS_H
#define ANISOTROPICELECTROSTATICS_H

#include "Kernel.h"

class AnisotropicElectrostatics;

template<>
InputParameters validParams<AnisotropicElectrostatics>();

class AnisotropicElectrostatics: public Kernel
{
public:

  AnisotropicElectrostatics(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const Real _inplane_permittivity;
  const Real _outofplane_permittivity;
  const Real _len_scale;

};
#endif
