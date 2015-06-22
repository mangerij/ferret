#ifndef PX_FIELDAUX_H
#define PX_FIELDAUX_H

#include "AuxKernel.h"


//Forward declarations
class Px_fieldAux;

template<>
InputParameters validParams<Px_fieldAux>();

class Px_fieldAux : public AuxKernel
{
public:
  Px_fieldAux( const std::string & name, InputParameters parameters );

  virtual ~Px_fieldAux() {}

protected:
  virtual Real computeValue();

private:
  const Real _permittivity_int;
  const Real _permittivity_ext;
 // const VariableGradient&  _potential_int_grad;
  const VariableGradient&  _potential_ext_grad;
};

#endif 


