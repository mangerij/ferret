/**
 * @file   EulerSkyrmionPCubeTerm.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 */

#ifndef EULERSKYRMIONPCUBETERM_H
#define EULERSKYRMIONPCUBETERM_H

#include "Kernel.h"

//Forward Declarations
class EulerSkyrmionPCubeTerm;

template<>
InputParameters validParams<EulerSkyrmionPCubeTerm>();

class EulerSkyrmionPCubeTerm: public Kernel
{
public:

  EulerSkyrmionPCubeTerm(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  const unsigned int _P_var;
  const VariableValue & _P;
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm
  const Real _P0;

};
#endif //EULERSKYRMIONPCUBETERM_H
