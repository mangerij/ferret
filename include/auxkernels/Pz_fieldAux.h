#ifndef PZ_FIELDAUX_H
#define PZ_FIELDAUX_H

#include "AuxKernel.h"


//Forward declarations
class Pz_fieldAux;

template<>
InputParameters validParams<Pz_fieldAux>();


class Pz_fieldAux : public AuxKernel
{
public:
  Pz_fieldAux( const std::string & name, InputParameters parameters );

  virtual ~Pz_fieldAux() {}

protected:
  virtual Real computeValue();

private:
  const Real _permittivity_int;
  const Real _permittivity_ext;
 // const VariableGradient&  _potential_int_grad;
  const VariableGradient&  _potential_ext_grad;

};

#endif 


