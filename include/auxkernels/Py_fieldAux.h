#ifndef PY_FIELDAUX_H
#define PY_FIELDAUX_H

#include "AuxKernel.h"


//Forward declarations
class Py_fieldAux;

template<>
InputParameters validParams<Py_fieldAux>();


class Py_fieldAux : public AuxKernel
{
public:
  Py_fieldAux( const std::string & name, InputParameters parameters );

  virtual ~Py_fieldAux() {}

protected:
  virtual Real computeValue();

private:
  const Real _permittivity_int;
  const Real _permittivity_ext;
 // const VariableGradient&  _potential_int_grad;
  const VariableGradient&  _potential_ext_grad;
};

#endif 


