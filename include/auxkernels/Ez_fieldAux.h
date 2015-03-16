#ifndef EZ_FIELDAUX_H
#define EZ_FIELDAUX_H

#include "AuxKernel.h"


//Forward declarations
class Ez_fieldAux;

template<>
InputParameters validParams<Ez_fieldAux>();


class Ez_fieldAux : public AuxKernel
{
public:
  Ez_fieldAux( const std::string & name, InputParameters parameters );

  virtual ~Ez_fieldAux() {}

protected:
  virtual Real computeValue();

private:
  const VariableGradient&  _potential_int_grad;
  const VariableGradient&  _potential_ext_grad;

};

#endif 


