/**
 * @file   WallEnergyDerivative.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu May 30 12:00:20 2013
 *
 * @brief
 *
 *
 */

#ifndef WALLENERGYDERIVATIVE_H
#define WALLENERGYDERIVATIVE_H

#include "Kernel.h"

class WallEnergyDerivative;

template<>
InputParameters validParams<WallEnergyDerivative>();

class WallEnergyDerivative: public Kernel
{
public:

  WallEnergyDerivative(const std::string & name, InputParameters parameters);

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
  const Real _G110,_G11, _G12, _G44, _G44P;
  const Real _len_scale;

};
#endif //WALLENERGYDERIVATIVE_H
