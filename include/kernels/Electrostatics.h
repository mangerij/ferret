/**
 * @file   Electrostatics.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jun 23 10:07:53 2015
 * @modified J. Mangeri <mangerij@anl.gov>
 * @brief  Laplacian operator with permittivity.
 *
 *
 */

#ifndef ELECTROSTATICS_H
#define ELECTROSTATICS_H

#include "Kernel.h"

class Electrostatics;

template<>
InputParameters validParams<Electrostatics>();

class Electrostatics: public Kernel
{
public:

  Electrostatics(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const Real _permittivity;
  const Real _len_scale;

};
#endif
