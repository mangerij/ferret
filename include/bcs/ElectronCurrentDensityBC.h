#ifndef ELECTRONCURRENTDENSITYBC_H
#define ELECTRONCURRENTDENSITYBC_H

#include "IntegratedBC.h"

class ElectronCurrentDensityBC : public IntegratedBC
{
public:
  ElectronCurrentDensityBC(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();

private:
const Real _Ec;
const Real _Nc;
const Real _T;
const Real _Kb;
const Real _q;
const Real _mun;
};

#endif
