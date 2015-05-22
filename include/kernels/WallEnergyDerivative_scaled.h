/**
 * @file   WallEnergyDerivative_scaled.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu May 30 12:00:20 2013
 *
 * @brief
 *
 *
 */

#ifndef WALLENERGYDERIVATIVE_SCALED_H
#define WALLENERGYDERIVATIVE_SCALED_H

#include "Kernel.h"

class WallEnergyDerivative_scaled;

template<>
InputParameters validParams<WallEnergyDerivative_scaled>();

class WallEnergyDerivative_scaled: public Kernel
{
public:

  WallEnergyDerivative_scaled(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableGradient& _polar_i_grad;
  const VariableGradient& _polar_j_grad;
  const VariableGradient& _polar_k_grad;
  const unsigned int _ii, _jj, _kk;
  const Real _alpha1;
  const Real _G110,_G11, _G12, _G44, _G44P;
  const Real _len_scale;
  const Real _energy_scale;

};
#endif //WALLENERGYDERIVATIVE_SCALED_H
