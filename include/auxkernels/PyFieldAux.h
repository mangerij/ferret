#ifndef PYFIELDAUX_H
#define PYFIELDAUX_H

#include "AuxKernel.h"


//Forward declarations
class PyFieldAux;

template<>
InputParameters validParams<PyFieldAux>();


class PyFieldAux : public AuxKernel
{
public:
  PyFieldAux( const std::string & name, InputParameters parameters );

  virtual ~PyFieldAux() {}

protected:
  virtual Real computeValue();

private:
  const Real _permittivity_int;
  const Real _permittivity_ext;
 // const VariableGradient&  _potential_int_grad;
  const VariableGradient & _potential_ext_grad;
};

#endif
