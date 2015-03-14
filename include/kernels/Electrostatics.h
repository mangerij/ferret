/**
 * @file   Electrostatics.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jun 11 10:07:53 2013
 *
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

  Electrostatics(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const Real _permittivity;
  const Real _len_scale;
  const Real _polar_electric_scale;

};
#endif
