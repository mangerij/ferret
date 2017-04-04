/**
 * @file   RenormalizedFreeEnergy.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Jun 1 12:00:20 2015
 *
 *
 */

#ifndef RENORMALIZEDFREEENERGY_H
#define RENORMALIZEDFREEENERGY_H

#include "Kernel.h"

class RenormalizedFreeEnergy;

template<>
InputParameters validParams<RenormalizedFreeEnergy>();

class RenormalizedFreeEnergy: public Kernel
{
public:

  RenormalizedFreeEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const Real _epsilon;
  const Real _T;
};
#endif //RENORMALIZEDFREEENERGY_H
