/****************************************************************/
/* Modified from MOOSE for FERRET module on Jun. 15 2015        */
/* to implement scaling on StressDivergenceTensorsScaled.C      */
/* @contact J. Mangeri <mangerij@anl.gov>                       */
/****************************************************************/

#ifndef TENSORMECHANICSACTIONSCALED_H
#define TENSORMECHANICSACTIONSCALED_H

#include "Action.h"

class TensorMechanicsActionScaled;

template<>
InputParameters validParams<TensorMechanicsActionScaled>();

class TensorMechanicsActionScaled : public Action
{
public:
  TensorMechanicsActionScaled(const InputParameters & parameters);

  virtual void act();

private:

};

#endif //TENSORMECHANICSACTIONSCALED_H
