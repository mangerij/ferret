


#ifndef BULKENERGYPSTO_H
#define BULKENERGYPSTO_H

//TODO: include the base header
#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class BulkEnergyPSTO;

template<>
InputParameters validParams<BulkEnergyPSTO>();

//TODO: change the base class!
class BulkEnergyPSTO : public ElementIntegralPostprocessor
{
public:
  BulkEnergyPSTO(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const Real _alpha1, _alpha2, _alpha3, _alpha4, _alpha5, _T, _Tc;

};

#endif
