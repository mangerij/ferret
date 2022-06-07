#ifndef CARRIERINT_H
#define CARRIERINT_H

#include "Kernel.h"

class CarrierInt: public Kernel
{
public:

  CarrierInt(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const Real _ni;

};
#endif //CARRIERINT_H
