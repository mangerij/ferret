#ifndef HOLES_H
#define HOLES_H

#include "AuxKernel.h"

class Holes: public AuxKernel
{
public:
  Holes(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~Holes() {}

protected:
  virtual Real computeValue();

private:
  const Real _Ev;
  const Real _Nv;
  const Real _T;
  const Real _Kb;
  const Real _q;
  const VariableValue & _potential_E_int;

};
#endif //HOLES_H
