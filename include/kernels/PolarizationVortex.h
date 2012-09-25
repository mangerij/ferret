#ifndef POLARIZATIONVORTEX_H
#define POLARIZATIONVORTEX_H

#include "Kernel.h"

//Forward Declarations
class PolarizationVortex;

template<>
InputParameters validParams<PolarizationVortex>();

class PolarizationVortex : public Kernel
{
public:

  PolarizationVortex(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  //virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _i;
  const std::string  _p;
  Real _a_x, _a_y, _c, _R,_L;
  VariableValue& _P_1, _P_2;
  bool _debug;
};
#endif //POLARIZATIONVORTEX_H
