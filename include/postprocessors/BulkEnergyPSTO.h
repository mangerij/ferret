/**
 * @file   BulkEnergy.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jun  4 15:05:42 2013
 *
 * @brief
 *
 *
 */


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
  const Real _alpha1, _alpha2, _alpha3, _alpha4, _alpha5,_x1, _x2, _x3, _x4, _x5, _x6, _epsilon, _T, _Tc;

};

#endif
