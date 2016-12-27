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
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const Real _alpha1, _alpha11, _alpha12, _alpha111, _alpha112,_alpha123;
  const Real _len_scale;
};

#endif // CUBICDIELECTRICTENSOR_H
