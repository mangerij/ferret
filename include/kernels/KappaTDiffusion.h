/**
 * @file   KappaTDiffusion.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * Note we hard-code the kappa T dependence here for PTO.
 *
 *
 */

#ifndef KAPPATDIFFUSION_H
#define KAPPATDIFFUSION_H

#include "Kernel.h"

class KappaTDiffusion;

template<>
InputParameters validParams<KappaTDiffusion>();

class KappaTDiffusion: public Kernel
{
public:

  KappaTDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const VariableValue & _temperature;
  const VariableGradient & _temperature_grad;
  const Real _c0, _c1, _c2, _c3, _c4;


};
#endif
