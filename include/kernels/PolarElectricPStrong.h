/**
 * @file   KernelTemplate.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu May 30 12:00:20 2013
 *
 * @brief
 *
 *
 */

#ifndef POLARELECTRICPSTRONG_H
#define POLARELECTRICPSTRONG_H

#include "Kernel.h"

class PolarElectricPStrong;

template<>
InputParameters validParams<PolarElectricPStrong>();

class PolarElectricPStrong: public Kernel
{
public:

  PolarElectricPStrong(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _component;
  const VariableGradient&  _potential_int_grad;
  const VariableGradient&  _potential_ext_grad;
  const Real _len_scale;
};
#endif //POLARELECTRICPSTRONG_H
