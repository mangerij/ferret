#ifndef MAGNETIZING_H
#define MAGNETIZING_H

#include "Kernel.h"

class Magnetizing;

template<>
InputParameters validParams<Magnetizing>();


class Magnetizing: public Kernel
{
public:

  Magnetizing(const std::string & name, InputParameters parameters);

protected:
  MaterialProperty<RealVectorValue >& _p;
  virtual Real computeQpResidual(); 
};

#endif //MAGNETIZING
