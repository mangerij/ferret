#ifndef POLARIZATIONVORTEXAUX_H
#define POLARIZATIONVORTEXAUX_H

#include "AuxKernel.h"
#include "FerretBase.h"

//Forward Declarations
class PolarizationVortexAux;

template<>
InputParameters validParams<PolarizationVortexAux>();

class PolarizationVortexAux : public FerretBase, public AuxKernel
{
public:

  PolarizationVortexAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

private:
  const unsigned int _i;
  const std::string  _p;
  Real _a_x, _a_y, _c, _R,_L;
  // VariableValue& _P_1, _P_2;
};
#endif //POLARIZATIONVORTEXAUX_H
