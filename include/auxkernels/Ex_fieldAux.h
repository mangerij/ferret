#ifndef EX_FIELDAUX_H
#define EX_FIELDAUX_H

#include "AuxKernel.h"


//Forward declarations
class Ex_fieldAux;

template<>
InputParameters validParams<Ex_fieldAux>();

class Ex_fieldAux : public AuxKernel
{
public:
  Ex_fieldAux( const std::string & name, InputParameters parameters );

  virtual ~Ex_fieldAux() {}

protected:
  virtual Real computeValue();

private:
  const VariableGradient&  _potential_int_grad;
  const VariableGradient&  _potential_ext_grad;
};

#endif 


