/**
 * @file   ElectrostrictiveEnergy.h
 * @author J. Mangeri <mangerij@anl.gov>
 *
 */

#ifndef ELECTROSTRICTIVEENERGY_H
#define ELECTROSTRICTIVEENERGY_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "ElementIntegralPostprocessor.h"
#include "ComputeElectrostrictiveTensor.h"
#include "ComputeEigenstrain.h"

//Forward Declarations
class ElectrostrictiveEnergy;

template<>
InputParameters validParams<ElectrostrictiveEnergy>();


class ElectrostrictiveEnergy : public ElementIntegralPostprocessor
{
public:
  ElectrostrictiveEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

private:
  const MaterialProperty<RankFourTensor> & _electrostrictive_tensor;
  const MaterialProperty<RankTwoTensor> & _eigenstrain;
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
