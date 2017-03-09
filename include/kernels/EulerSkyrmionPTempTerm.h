/**
 * @file   EulerSkyrmionPTempTerm.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 */

#ifndef EULERSKYRMIONPTEMPTERM_H
#define EULERSKYRMIONPTEMPTERM_H

#include "Kernel.h"

//Forward Declarations
class EulerSkyrmionPTempTerm;

template<>
InputParameters validParams<EulerSkyrmionPTempTerm>();

class EulerSkyrmionPTempTerm: public Kernel
{
public:

  EulerSkyrmionPTempTerm(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _theta_var;
  const unsigned int _P_var;
  const VariableValue & _theta;
  const VariableGradient & _theta_grad;
  const VariableValue & _P;
  const Real _t;
  const Real _kappa;
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm
  const Real _xi0;

};
#endif //EULERSKYRMIONPTEMPTERM_H
