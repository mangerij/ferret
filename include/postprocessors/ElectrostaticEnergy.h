/**
 * @file   ElectrostaticEnergy.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jul 30 17:08:18 2013
 *
 * @brief compute the electric Energy 0.5*D\dot E
 *
 *
 */

#ifndef ELECTROSTATICENERGY_H
#define ELECTROSTATICENERGY_H

//TODO: include the base header
#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class ElectrostaticEnergy;

template<>
InputParameters validParams<ElectrostaticEnergy>();

//TODO: change the base class!
class ElectrostaticEnergy : public ElementIntegralPostprocessor
{
public:
  ElectrostaticEnergy(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const VariableGradient& _potential_grad;   //total potential
  const VariableGradient& _potential_ext_grad;   //for external electric field
  const Real _permittivity;
  const Real _len_scale;
};

#endif
