/**
 * @file   CubicDielectricTensor.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief calculate the components of the anisotropic dielectric tensor
 *        assuming cubic symmetry of the parent phase
 */

#ifndef CUBICDIELECTRICTENSOR_H
#define CUBICDIELECTRICTENSOR_H

#include "AuxKernel.h"
#include "ComputeElectrostrictiveTensor.h"

//Forward Declarations
class CubicDielectricTensor;

template<>
InputParameters validParams<CubicDielectricTensor>();

class CubicDielectricTensor : public AuxKernel
{
public:
  CubicDielectricTensor(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  const MaterialProperty<RankFourTensor> & _electrostrictive_tensor;
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const VariableGradient & _disp_x_grad;
  const VariableGradient & _disp_y_grad;
  const VariableGradient & _disp_z_grad;
  const Real _alpha1, _alpha11, _alpha12, _alpha111, _alpha112,_alpha123;
  const unsigned int _first_deriv, _second_deriv;
  const Real _len_scale;
};

#endif // CUBICDIELECTRICTENSOR_H
