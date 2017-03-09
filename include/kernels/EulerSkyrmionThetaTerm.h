/**
 * @file   EulerSkyrmionThetaTerm.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 */

#ifndef EULERSKYRMIONTHETATERM_H
#define EULERSKYRMIONTHETATERM_H

#include "Kernel.h"

//Forward Declarations
class EulerSkyrmionThetaTerm;

template<>
InputParameters validParams<EulerSkyrmionThetaTerm>();

class EulerSkyrmionThetaTerm: public Kernel
{
public:

  EulerSkyrmionThetaTerm(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _theta_var;
  const unsigned int _P_var;
  const VariableValue & _theta;
  const VariableSecond & _second_u;
  const VariableTestSecond & _second_test;
  const VariablePhiSecond & _second_phi;
  const VariableValue & _P;
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm
  const Real _xi0;

};
#endif //EULERSKYRMIONTHETATERM_H
