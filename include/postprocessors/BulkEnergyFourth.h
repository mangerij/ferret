/**
 * @file   BulkEnergyFourth.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Aug  2 15:05:42 2015
 *
 * @brief
 *
 *
 */


#ifndef BULKENERGYFOURTH_H
#define BULKENERGYFOURTH_H

//TODO: include the base header
#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class BulkEnergyFourth;

template<>
InputParameters validParams<BulkEnergyFourth>();

//TODO: change the base class!
class BulkEnergyFourth : public ElementIntegralPostprocessor
{
public:
  BulkEnergyFourth(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const Real _alpha1, _alpha11, _alpha12;
  const Real _len_scale;
};

#endif
