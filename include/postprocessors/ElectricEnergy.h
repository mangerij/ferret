/**
 * @file   ElectricEnergy.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jul 30 17:08:18 2013
 *
 * @brief compute the electric Energy 0.5*D\dot E
 *
 *
 */

#ifndef ELECTRICENERGY_H
#define ELECTRICENERGY_H

//TODO: include the base header
#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class ElectricEnergy;

template<>
InputParameters validParams<ElectricEnergy>();

//TODO: change the base class!
class ElectricEnergy : public ElementIntegralPostprocessor
{
public:
  ElectricEnergy(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const VariableGradient& _potential_grad;
  const Real _permittivity;
  const Real _len_scale;
};

#endif
