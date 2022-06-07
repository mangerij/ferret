#ifndef HOLECURRENTDENSITYBC_H
#define HOLECURRENTDENSITYBC_H

#include "IntegratedBC.h"

class HoleCurrentDensityBC : public IntegratedBC
{
public:
  HoleCurrentDensityBC(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();

private:
const Real _Ev;
const Real _Nv;
const Real _T;
const Real _Kb;
const Real _q;
const Real _mup;
};

#endif
