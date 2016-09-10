#ifndef MIEELECFIELDREALS_H
#define MIEELECFIELDREALS_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class MieElecFieldReals;

template<>
InputParameters validParams<MieElecFieldReals>();


class MieElecFieldReals : public AuxKernel
{
public:
  MieElecFieldReals(const InputParameters & parameters);

  virtual ~MieElecFieldReals() {}

protected:
  virtual Real computeValue();
  const Real _a, _omega, _c, _epsilonI, _sigmaI, _epsilonII, _sigmaII, _L, _nh, _order, _scale, _component;
private:

};

#endif // MIEELECFIELDREALS_H
