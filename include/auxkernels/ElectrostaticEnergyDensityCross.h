/**
 * @file   ElectrostaticEnergyDensityCross.h   ElectrostaticEnergyDensityCross.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Oct  2 17:13:46 2013   Wed Oct  2 17:12:30 2013
 * @brief  ElectroEnergyDensity Cross Term
 */

#ifndef ELECTROSTATICENERGYDENSITYCROSS_H
#define ELECTROSTATICENERGYDENSITYCROSS_H

#include "AuxKernel.h"

//Forward Declarations
class ElectrostaticEnergyDensityCross;

template<>
InputParameters validParams<ElectrostaticEnergyDensityCross>();

class ElectrostaticEnergyDensityCross : public AuxKernel
{
public:
  ElectrostaticEnergyDensityCross(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();
  const VariableGradient& _potential_grad;
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const Real _len_scale;
  const Real _energy_scale;
};

#endif // ELECTROSTATICENERGYDENSITYCROSS_H
