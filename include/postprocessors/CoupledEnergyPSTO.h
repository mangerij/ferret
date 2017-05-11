


#ifndef COUPLEDENERGYPSTO_H
#define COUPLEDENERGYPSTO_H

//TODO: include the base header
#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class CoupledEnergyPSTO;

template<>
InputParameters validParams<CoupledEnergyPSTO>();

//TODO: change the base class!
class CoupledEnergyPSTO : public ElementIntegralPostprocessor
{
public:
  CoupledEnergyPSTO(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const Real _x1, _x2, _x3, _x4, _x5, _x6, _epsilon;

};

#endif
