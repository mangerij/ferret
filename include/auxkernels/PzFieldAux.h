#ifndef PZFIELDAUX_H
#define PZFIELDAUX_H

#include "AuxKernel.h"


//Forward declarations
class PzFieldAux;

template<>
InputParameters validParams<PzFieldAux>();


class PzFieldAux : public AuxKernel
{
public:
  PzFieldAux( const std::string & name, InputParameters parameters );

  virtual ~PzFieldAux() {}

protected:
  virtual Real computeValue();

private:
  const Real _permittivity_int;
  const Real _permittivity_ext;
 // const VariableGradient&  _potential_int_grad;
  const VariableGradient & _potential_ext_grad;
};

#endif
