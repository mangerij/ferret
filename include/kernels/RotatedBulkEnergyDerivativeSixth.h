/**
 * @file   RotatedBulkEnergyDerivativeSixth.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @date   Jun 1 12:00:20 2015
 *
 * @brief
 * This rotates the bulk energy functional. This should be the default Kernel to use but might be slower due to extra steps
 *
 */

#ifndef ROTATEDBULKENERGYDERIVATIVESIXTH_H
#define ROTATEDBULKENERGYDERIVATIVESIXTH_H

#include "Kernel.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"

class RotatedBulkEnergyDerivativeSixth;

template<>
InputParameters validParams<RotatedBulkEnergyDerivativeSixth>();

class RotatedBulkEnergyDerivativeSixth: public Kernel
{
public:

  RotatedBulkEnergyDerivativeSixth(const InputParameters & parameters);

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
  RealVectorValue _Euler_angles;
  const Real _alpha1, _alpha11, _alpha12, _alpha111, _alpha112,_alpha123;
  const Real _len_scale;
};
#endif //ROTATEDBULKENERGYDERIVATIVESIXTH_H
