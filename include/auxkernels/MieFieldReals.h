#ifndef MIEFIELDREALS_H
#define MIEFIELDREALS_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class MieFieldReals;

template<>
InputParameters validParams<MieFieldReals>();


class MieFieldReals : public AuxKernel
{
public:
  MieFieldReals(const InputParameters & parameters);

  virtual ~MieFieldReals() {}

protected:
  virtual Real computeValue();
  const Real _a, _omega, _c, _epsilonI, _sigmaI, _epsilonII, _sigmaII, _L, _order, _component;
private:

};

#endif // MIEFIELDREALS_H
