#ifndef ELECCURRENT_H
#define ELECCURRENT_H

#include "Kernel.h"

class ElecCurrent: public Kernel
{
public:

  ElecCurrent(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const Real _mun;
  const Real _Nc;
  const Real _q;
  const Real _Ec;
  const Real _Kb;
  const Real _T;
};
#endif //ELECCURRENT_H
