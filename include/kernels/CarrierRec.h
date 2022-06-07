#ifndef CARRIERREC_H
#define CARRIERREC_H

#include "Kernel.h"

class CarrierRec: public Kernel
{
public:

  CarrierRec(const InputParameters & parameters);

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
  const Real _b;

};
#endif //CARRIERCON_H
