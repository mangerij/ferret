/**
 * @file   BulkEnergyDerivativeSixth.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu May 30 12:00:20 2013
 *
 * @brief
 *
 *
 */

#ifndef BULKENERGYDERIVATIVESIXTH_H
#define BULKENERGYDERIVATIVESIXTH_H

#include "Kernel.h"

class BulkEnergyDerivativeSixth;

template<>
InputParameters validParams<BulkEnergyDerivativeSixth>();

class BulkEnergyDerivativeSixth: public Kernel
{
public:

  BulkEnergyDerivativeSixth(const InputParameters & parameters);

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
  const Real _alpha1, _alpha11, _alpha12, _alpha111, _alpha112,_alpha123, _kappa;
  const Real _len_scale;
};
#endif //BULKENERGYDERIVATIVESIXTH_H
