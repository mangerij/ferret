/**
 * @file   DielectricTensor.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief calculate the components of the anisotropic dielectric tensor
 */

#ifndef DIELECTRICTENSOR_H
#define DIELECTRICTENSOR_H

#include "AuxKernel.h"
#include "ComputeElectrostrictiveTensor.h"

//Forward Declarations
class DielectricTensor;
class RankTwoTensor;

template<>
InputParameters validParams<DielectricTensor>();

class DielectricTensor : public AuxKernel
{
public:
  DielectricTensor(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  const MaterialProperty<RankFourTensor> & _electrostrictive_tensor;
  const MaterialProperty<RankTwoTensor> & _elastic_strain;
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

#endif // DIELECTRICTENSOR_H
