#ifndef LC_H
#define LC_H

#include "Kernel.h"

// #include "Material.h"
// #include "DerivativeMaterialInterface.h"

//Forward Declarations
class Lc;

template<>
InputParameters validParams<Lc>();

class Lc: public Kernel
{
public:

  Lc(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const unsigned int _c_var;
  const VariableValue & _c;
  const Real _kappa;
};
#endif //Lc_H
