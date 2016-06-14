/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef MODIFIEDSTRESSDIVERGENCETENSORS_H
#define MODIFIEDSTRESSDIVERGENCETENSORS_H

#include "Kernel.h"
#include "RankFourTensor.h"
#include "ComputeElectrostrictiveTensor.h"
#include "RankTwoTensor.h"

//Forward Declarations
class ModifiedStressDivergenceTensors;
class RankFourTensor;
class RankTwoTensor;

template<>
InputParameters validParams<ModifiedStressDivergenceTensors>();

/**
 * StressDivergenceTensors mostly copies from StressDivergence.  There are small changes to use
 * RankFourTensor and RankTwoTensors instead of SymmElasticityTensors and SymmTensors.  This is done
 * to allow for more mathematical transparancy.
 */
class ModifiedStressDivergenceTensors : public Kernel
{
public:
  ModifiedStressDivergenceTensors(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  std::string _base_name;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;
  // MaterialProperty<RankTwoTensor> & _d_stress_dT;

  const unsigned int _component;

  /// Coupled displacement variables
  unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

  const bool _temp_coupled;

  const unsigned int _temp_var;

private:
  const MaterialProperty<RankFourTensor> & _electrostrictive_tensor;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const VariableGradient & _polar_x_grad;
  const VariableGradient & _polar_y_grad;
  const VariableGradient & _polar_z_grad;
  const Real _len_scale;

};

#endif //MODIFIEDSTRESSDIVERGENCETENSORS_H
