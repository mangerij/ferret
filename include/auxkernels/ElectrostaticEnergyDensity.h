/**
 * @file   ElectrostaticEnergyDensity.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Oct  2 17:13:46 2013
 * @brief  ElectroEnergyDensity
 */

#ifndef ELECTROSTATICENERGYDENSITY_H
#define ELECTROSTATICENERGYDENSITY_H

#include "AuxKernel.h"

//Forward Declarations
class ElectrostaticEnergyDensity;

template<>
InputParameters validParams<ElectrostaticEnergyDensity>();

class ElectrostaticEnergyDensity : public AuxKernel
{
public:
  ElectrostaticEnergyDensity(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();
  const VariableGradient& _potential_int_grad;
  const VariableGradient& _potential_ext_grad;
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const Real _len_scale;
  const Real _energy_scale;
};

#endif // ELECTROSTATICENERGYDENSITY_H
