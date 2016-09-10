#ifndef MIEELECFIELDIMAG_H
#define MIEELECFIELDIMAG_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class MieElecFieldImag;

template<>
InputParameters validParams<MieElecFieldImag>();


class MieElecFieldImag : public AuxKernel
{
public:
  MieElecFieldImag(const InputParameters & parameters);

  virtual ~MieElecFieldImag() {}

protected:
  virtual Real computeValue();
  const Real _a, _omega, _c, _epsilonI, _sigmaI, _epsilonII, _sigmaII, _L, _nh, _order, _scale, _component;
private:

};

#endif // MIEELECFIELDIMAG_H
