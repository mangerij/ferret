#ifndef CARRIERCON_H
#define CARRIERCON_H

#include "Kernel.h"

class CarrierCon: public Kernel
{
public:

  CarrierCon(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const Real _Ev;
  const Real _Ec;
  const Real _Nv;
  const Real _Nc;
  const Real _T;
  const Real _Kb;
  const Real _q;
  const Real _Na;
  const Real _Nd;
};
#endif //CARRIERCON_H
