/************************************************************************/
/* ScreenedBC:                                                          */
/*     This BC is intended to implement a \lambda P                     */
/*     + \epsilon * _grad\phi = 0 at a boundary where                   */        
/*     \lambda is a screening param.                                    */
/************************************************************************/


#ifndef SCREENEDBC_H
#define SCREENEDBC_H

#include "IntegratedBC.h"

//Forward Declarations
class ScreenedBC;

template<>
InputParameters validParams<ScreenedBC>();

class ScreenedBC : public IntegratedBC
{
public:

  ScreenedBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

private:
  const unsigned int _component;
  const VariableGradient &  _potential_int_grad;
  const Real _permittivity;
  const Real _lambda;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
};

#endif /* SCREENEDBC_H */
