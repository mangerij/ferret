/**
 * @file   ElectricStatics.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jun 11 10:07:53 2013
 * 
 * @brief  Laplacian operator with permittivity.
 * 
 * 
 */

#ifndef ELECTRICSTATICS_H
#define ELECTRICSTATICS_H

#include "Kernel.h"

class ElectricStatics;

template<>
InputParameters validParams<ElectricStatics>();

class ElectricStatics: public Kernel
{
public:

  ElectricStatics(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

private:
  
  const Real _permittivity;
  
};
#endif //ELECTRICDTSTATICS_H
