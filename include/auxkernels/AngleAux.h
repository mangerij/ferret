#ifndef ANGLEAUX_H
#define ANGLEAUX_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

//Forward declarations
class AngleAux;

template<>
InputParameters validParams<AngleAux>();


class AngleAux : public AuxKernel
{
public:
  AngleAux(const InputParameters & parameters);

  virtual ~AngleAux() {}

protected:
  virtual Real computeValue();
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
};

#endif // ANGLEAUX_H
