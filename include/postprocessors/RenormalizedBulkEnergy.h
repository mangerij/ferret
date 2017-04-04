/**
 * @file   RenormalizedBulkEnergy.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief
 *
 *
 */


#ifndef RENORMALIZEDBULKENERGY_H
#define RENORMALIZEDBULKENERGY_H

//TODO: include the base header
#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class RenormalizedBulkEnergy;

template<>
InputParameters validParams<RenormalizedBulkEnergy>();

//TODO: change the base class!
class RenormalizedBulkEnergy : public ElementIntegralPostprocessor
{
public:
  RenormalizedBulkEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const Real _epsilon;
  const Real _T;

};

#endif
