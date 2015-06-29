#ifndef PXFIELDAUX_H
#define PXFIELDAUX_H

#include "AuxKernel.h"


//Forward declarations
class PxFieldAux;

template<>
InputParameters validParams<PxFieldAux>();

class PxFieldAux : public AuxKernel
{
public:
  PxFieldAux( const std::string & name, InputParameters parameters );

  virtual ~PxFieldAux() {}

protected:
  virtual Real computeValue();

private:
  const Real _permittivity_int;
  const Real _permittivity_ext;
 // const VariableGradient&  _potential_int_grad;
  const VariableGradient & _potential_ext_grad;
};

#endif
