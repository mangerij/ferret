#ifndef ELECTRONS_H
#define ELECTRONS_H

#include "AuxKernel.h"

class Electrons: public AuxKernel
{
public:
  Electrons(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~Electrons() {}

protected:
  virtual Real computeValue();

private:
  const Real _Ec;
  const Real _Nc;
  const Real _T;
  const Real _Kb;
  const Real _q;
  const VariableValue & _potential_E_int;

};
#endif //ELECTRONS_H
