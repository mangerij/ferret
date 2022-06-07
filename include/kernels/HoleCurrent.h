#ifndef HOLECURRENT_H
#define HOLECURRENT_H

#include "Kernel.h"

class HoleCurrent: public Kernel
{
public:

  HoleCurrent(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const Real _Ev;
  const Real _Nv;
  const Real _T;
  const Real _Kb;
  const Real _q;
  const Real _mup;

};
#endif //HOLECURRENT_H
