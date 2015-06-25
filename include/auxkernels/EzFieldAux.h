#ifndef EZFIELDAUX_H
#define EZFIELDAUX_H

#include "AuxKernel.h"


//Forward declarations
class EzFieldAux;

template<>
InputParameters validParams<EzFieldAux>();


class EzFieldAux : public AuxKernel
{
public:
  EzFieldAux( const std::string & name, InputParameters parameters );

  virtual ~EzFieldAux() {}

protected:
  virtual Real computeValue();

private:
  const VariableGradient & _potential_int_grad;
  const VariableGradient & _potential_ext_grad;

};

#endif /* EZFIELDAUX_H */
