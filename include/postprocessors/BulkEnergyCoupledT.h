/**
 * @file   BulkEnergyCoupledT.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 * @brief
 *
 *
 */


#ifndef BULKENERGYCOUPLEDT_H
#define BULKENERGYCOUPLEDT_H

#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class BulkEnergyCoupledT;

template<>
InputParameters validParams<BulkEnergyCoupledT>();

class BulkEnergyCoupledT : public ElementIntegralPostprocessor
{
public:
  BulkEnergyCoupledT(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const VariableValue & _temperature;
  const Real _alpha0, _alpha11, _alpha12, _alpha111, _alpha112, _alpha123, _Tc;
  const Real _len_scale;

};

#endif
