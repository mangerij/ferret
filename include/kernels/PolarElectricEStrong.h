/**
 * @file   KernelTemplate.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu May 30 12:00:20 2013
 *
 * @brief
 *
 *
 */

#ifndef POLARELECTRICESTRONG_H
#define POLARELECTRICESTRONG_H

#include "Kernel.h"

class PolarElectricEStrong;

template<>
InputParameters validParams<PolarElectricEStrong>();

class PolarElectricEStrong: public Kernel
{
public:

  PolarElectricEStrong(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
   const Real _permittivity;
   const unsigned int _polar_x_var;
   const unsigned int _polar_y_var;
   const unsigned int _polar_z_var;
   const VariableValue & _polar_x;
   const VariableValue & _polar_y;
   const VariableValue & _polar_z;
   const Real _polar_scale;
   const Real _len_scale;

};
#endif //POLARELECTRICESTRONG_H
