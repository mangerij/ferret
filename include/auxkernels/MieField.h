#ifndef MIEFIELD_H
#define MIEFIELD_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class MieField;

template<>
InputParameters validParams<MieField>();


class MieField : public AuxKernel
{
public:
  MieField(const InputParameters & parameters);

  virtual ~MieField() {}

protected:
  virtual Real computeValue();
  const Real _a,_L, _order, _component;
private:

};

#endif // MIEFIELD_H
