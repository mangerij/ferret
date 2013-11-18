/**
 * @file   BulkEnergyDensity.h   BulkEnergyDensity.h
 * @author S. Gu <sgu@anl.gov>
 * @brief  calculate the bulk energy density:
 *
 *
 */

#ifndef BULKENERGYDENSITY_H
#define BULKENERGYDENSITY_H

#include "AuxKernel.h"


//Forward Declarations
class BulkEnergyDensity;

template<>
InputParameters validParams<BulkEnergyDensity>();

class BulkEnergyDensity : public AuxKernel
{
public:
  BulkEnergyDensity(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const Real _alpha1, _alpha11, _alpha12, _alpha111, _alpha112,_alpha123;
  const Real _len_scale;
  const Real _energy_scale;
};

#endif // BULKENERGYDENSITY_H
