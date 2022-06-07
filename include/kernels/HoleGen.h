#ifndef HOLEGEN_H
#define HOLEGEN_H

#include "Kernel.h"

class HoleGen: public Kernel
{
public:

  HoleGen(const InputParameters & parameters);

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
#endif //HOLEGEN_H
