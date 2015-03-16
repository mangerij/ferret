#ifndef EY_FIELDAUX_H
#define EY_FIELDAUX_H

#include "AuxKernel.h"


//Forward declarations
class Ey_fieldAux;

template<>
InputParameters validParams<Ey_fieldAux>();


class Ey_fieldAux : public AuxKernel
{
public:
  Ey_fieldAux( const std::string & name, InputParameters parameters );

  virtual ~Ey_fieldAux() {}

protected:
  virtual Real computeValue();

private:
  const VariableGradient&  _potential_int_grad;
  const VariableGradient&  _potential_ext_grad;
};

#endif 


