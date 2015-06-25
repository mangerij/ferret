#ifndef EXFIELDAUX_H
#define EXFIELDAUX_H

#include "AuxKernel.h"


//Forward declarations
class ExFieldAux;

template<>
InputParameters validParams<ExFieldAux>();

class ExFieldAux : public AuxKernel
{
public:
  ExFieldAux( const std::string & name, InputParameters parameters );

  virtual ~ExFieldAux() {}

protected:
  virtual Real computeValue();

private:
  const VariableGradient & _potential_int_grad;
  const VariableGradient & _potential_ext_grad;
};

#endif /* EXFIELDAUX_H */
