/**
 * @file   BulkEnergy.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jun  4 15:05:42 2013
 * 
 * @brief  
 * 
 * 
 */


#ifndef BULKENERGY_H
#define BULKENERGY_H

//TODO: include the base header
#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class BulkEnergy;

template<>
InputParameters validParams<BulkEnergy>();

//TODO: change the base class!
class BulkEnergy : public ElementIntegralPostprocessor  
{
public:
  BulkEnergy(const std::string & name, InputParameters parameters);
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const Real _alpha1, _alpha11, _alpha12, _alpha111, _alpha112,_alpha123; 

protected:
  virtual Real computeQpIntegral();
};

#endif
