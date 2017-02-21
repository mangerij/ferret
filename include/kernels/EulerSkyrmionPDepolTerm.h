/**
 * @file   EulerSkyrmionPDepolTerm.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 */

#ifndef EULERSKYRMIONPDEPOLTERM_H
#define EULERSKYRMIONPDEPOLTERM_H

#include "Kernel.h"

//Forward Declarations
class EulerSkyrmionPDepolTerm;

template<>
InputParameters validParams<EulerSkyrmionPDepolTerm>();

class EulerSkyrmionPDepolTerm: public Kernel
{
public:

  EulerSkyrmionPDepolTerm(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
private:
  const unsigned int _theta_var;
  const VariableValue & _theta;
  const Real _edep;
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm


};
#endif //EULERSKYRMIONPDEPOLTERM_H
