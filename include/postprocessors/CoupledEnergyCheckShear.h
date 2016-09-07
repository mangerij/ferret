/**
 * @file   CoupledEnergyCheckShear.h
 */

#ifndef COUPLEDENERGYCHECKSHEAR_H
#define COUPLEDENERGYCHECKSHEAR_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "ElementIntegralPostprocessor.h"
#include "ComputeElectrostrictiveTensor.h"
#include "ComputeEigenstrain.h"

//Forward Declarations
class CoupledEnergyCheckShear;

template<>
InputParameters validParams<CoupledEnergyCheckShear>();


class CoupledEnergyCheckShear : public ElementIntegralPostprocessor
{
public:
  CoupledEnergyCheckShear(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

private:
  const MaterialProperty<RankFourTensor> & _electrostrictive_tensor;
  const MaterialProperty<RankTwoTensor> & _stress_free_strain;
  const VariableGradient & _disp_x_grad;
  const VariableGradient & _disp_y_grad;
  const VariableGradient & _disp_z_grad;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const Real _artificial;
  const Real _len_scale;
};

#endif
