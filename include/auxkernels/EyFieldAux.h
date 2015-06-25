#ifndef EYFIELDAUX_H
#define EYFIELDAUX_H

#include "AuxKernel.h"


//Forward declarations
class EyFieldAux;

template<>
InputParameters validParams<EyFieldAux>();


class EyFieldAux : public AuxKernel
{
public:
  EyFieldAux( const std::string & name, InputParameters parameters );

  virtual ~EyFieldAux() {}

protected:
  virtual Real computeValue();

private:
  const VariableGradient & _potential_int_grad;
  const VariableGradient & _potential_ext_grad;
};

#endif /* EYFIELDAUX_H */
