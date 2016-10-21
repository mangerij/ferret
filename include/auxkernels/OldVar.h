#ifndef OLDVAR_H
#define OLDVAR_H

#include "AuxKernel.h"

//Forward declarations
class OldVar;

template<>
InputParameters validParams<OldVar>();


class OldVar : public AuxKernel
{
public:
  OldVar(const InputParameters & parameters);

  virtual ~OldVar() {}

protected:
  virtual Real computeValue();

};

#endif // BANDGAP_H
