/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
/* Modified(Scaled) for FERRET, added time_scale  */

#ifndef TIMEDERIVATIVE_SCALED_H
#define TIMEDERIVATIVE_SCALED_H

#include "TimeKernel.h"

// Forward Declaration
class TimeDerivativeScaled;

template<>
InputParameters validParams<TimeDerivativeScaled>();

class TimeDerivativeScaled : public TimeKernel
{
public:
  TimeDerivativeScaled(const std::string & name, InputParameters parameters);

  virtual void computeJacobian();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  bool _lumping;
  const Real _time_scale;
};

#endif //TIMEDERIVATIVE_SCALED_H
