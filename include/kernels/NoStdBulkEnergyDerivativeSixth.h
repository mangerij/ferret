/**
 * @file   NoStdBulkEnergyDerivativeSixth.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Jun 1 12:00:20 2015
 *
 * @brief
 *
 *
 */

#ifndef NOSTDBULKENERGYDERIVATIVESIXTH_H
#define NOSTDBULKENERGYDERIVATIVESIXTH_H

#include "Kernel.h"

class NoStdBulkEnergyDerivativeSixth;

template<>
InputParameters validParams<NoStdBulkEnergyDerivativeSixth>();

class NoStdBulkEnergyDerivativeSixth: public Kernel
{
public:

  NoStdBulkEnergyDerivativeSixth(const InputParameters & parameters);

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
  const Real _alpha1, _alpha11, _alpha12, _alpha111, _alpha112,_alpha123;
  const Real _len_scale;
};
#endif //NOSTDBULKENERGYDERIVATIVESIXTH_H
