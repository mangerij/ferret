/**
 * @file   PZTWallEnergyDerivative.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Dec 12 11:59:56 2015
 *
 * @brief
 *
 *
 */

#ifndef PZTWALLENERGYDERIVATIVE_H
#define PZTWALLENERGYDERIVATIVE_H

#include "Kernel.h"

class PZTWallEnergyDerivative;

template<>
InputParameters validParams<PZTWallEnergyDerivative>();

class PZTWallEnergyDerivative: public Kernel
{
public:

  PZTWallEnergyDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableGradient & _polar_i_grad;
  const VariableGradient & _polar_j_grad;
  const VariableGradient & _polar_k_grad;
  const unsigned int _ii, _jj, _kk;
  const Real _G110;
  const Real _len_scale;

};
#endif //PZTWALLENERGYDERIVATIVE_H
