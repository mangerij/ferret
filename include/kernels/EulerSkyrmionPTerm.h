/**
 * @file   EulerSkyrmionPTerm.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 */

#ifndef EULERSKYRMIONPTERM_H
#define EULERSKYRMIONPTERM_H

#include "Kernel.h"

//Forward Declarations
class EulerSkyrmionPTerm;

template<>
InputParameters validParams<EulerSkyrmionPTerm>();

class EulerSkyrmionPTerm: public Kernel
{
public:

  EulerSkyrmionPTerm(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

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
#endif //EULERSKYRMIONPTERM_H
