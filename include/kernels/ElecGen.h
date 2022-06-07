#ifndef ELECGEN_H
#define ELECGEN_H

#include "Kernel.h"

class ElecGen: public Kernel
{
public:

  ElecGen(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const Real _Ec;
  const Real _Nc;
  const Real _T;
  const Real _Kb;
  const Real _q;
  const Real _mun;

};
#endif //ELECGEN_H
