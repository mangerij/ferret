/**
 * @file   ThermalEnergy.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 */

#ifndef THERMALENERGY_H
#define THERMALENERGY_H

#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class ThermalEnergy;

template<>
InputParameters validParams<ThermalEnergy>();

class ThermalEnergy : public ElementIntegralPostprocessor
{
public:
  ThermalEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue & _temperature;
  const Real _const;
};

#endif
