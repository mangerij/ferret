/**
 * @file   CustomDT.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Fri Nov 15 13:17:57 2013
 *
 * @brief  Demonstrate the "ungluing of linear and nonlinear residual
 *
 *
 */
#ifndef CUSTOMDT_H_
#define CUSTOMDT_H_

#include "TimeStepper.h"
#include "libmesh/numeric_vector.h"
class CustomDT;

template<>
InputParameters validParams<CustomDT>();

class CustomDT :
  public TimeStepper
{
public:
  CustomDT(const std::string & name, InputParameters parameters);
protected:
  virtual Real computeInitialDT();
  virtual Real computeDT();
  Real _increase_rate;
};

#endif /* CustomDT_H_ */
