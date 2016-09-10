#ifndef MIEMAGFIELDREALS_H
#define MIEMAGFIELDREALS_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class MieMagFieldReals;

template<>
InputParameters validParams<MieMagFieldReals>();


class MieMagFieldReals : public AuxKernel
{
public:
  MieMagFieldReals(const InputParameters & parameters);

  virtual ~MieMagFieldReals() {}

protected:
  virtual Real computeValue();
  const Real _a, _omega, _c, _epsilonI, _sigmaI, _epsilonII, _sigmaII, _L, _nh, _order, _component;
private:

};

#endif // MIEMAGFIELDREALS_H
