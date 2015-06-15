/****************************************************************/
/* Modified from MOOSE for FERRET module on Jun. 15 2015        */
/* to implement scaling on StressDivergenceTensorsScaled.C      */
/* @contact J. Mangeri <mangerij@anl.gov>                       */
/****************************************************************/

#ifndef TENSORMECHANICSACTION_SCALED_H
#define TENSORMECHANICSACTION_SCALED_H

#include "Action.h"

class TensorMechanicsAction_scaled;

template<>
InputParameters validParams<TensorMechanicsAction_scaled>();

class TensorMechanicsAction_scaled : public Action
{
public:
  TensorMechanicsAction_scaled(const std::string & name, InputParameters params);

  virtual void act();

private:

};

#endif //TENSORMECHANICSACTION_H
