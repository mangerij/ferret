/**
 * @file   EulerSkyrmionThetaDepolTerm.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 */

#ifndef EULERSKYRMIONTHETADEPOLTERM_H
#define EULERSKYRMIONTHETADEPOLTERM_H

#include "Kernel.h"

//Forward Declarations
class EulerSkyrmionThetaDepolTerm;

template<>
InputParameters validParams<EulerSkyrmionThetaDepolTerm>();

class EulerSkyrmionThetaDepolTerm: public Kernel
{
public:

  EulerSkyrmionThetaDepolTerm(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const unsigned int _theta_var;
  const VariableValue & _theta;
  const Real _edep;
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm


};
#endif //EULERSKYRMIONTHETADEPOLTERM_H
