/**
 * @file   WallEnergyDerivative.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * Rotates WallEnergyDerivative
 *
 */

#ifndef ROTATEDWALLENERGYDERIVATIVE_H
#define ROTATEDWALLENERGYDERIVATIVE_H

#include "Kernel.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"

class RotatedWallEnergyDerivative;

template<>
InputParameters validParams<RotatedWallEnergyDerivative>();

class RotatedWallEnergyDerivative: public Kernel
{
public:

  RotatedWallEnergyDerivative(const InputParameters & parameters);

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
  RealVectorValue _Euler_angles;
  const unsigned int _ii, _jj, _kk;
  const Real _G110, _G11, _G12, _G44, _G44P;
  const Real _len_scale;

};
#endif //ROTATEDWALLENERGYDERIVATIVE_H
