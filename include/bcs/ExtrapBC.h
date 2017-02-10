/************************************************************************/
/* ExtrapBC:                                                            */
/*     This BC is intended to implement P + \beta dP/dn = 0             */
/*                                                                      */
/************************************************************************/



#ifndef EXTRAPBC_H
#define EXTRAPBC_H

#include "IntegratedBC.h"

//Forward Declarations
class ExtrapBC;

template<>
InputParameters validParams<ExtrapBC>();

class ExtrapBC : public IntegratedBC
{
public:

  ExtrapBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

private:
  const Real _beta;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const VariableGradient & _polar_z_grad;

};

#endif /* EXTRAPBC_H */
