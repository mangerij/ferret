#ifndef POLARIZATIONVORTEXACTION_H
#define POLARIZATIONVORTEXACTION_H

#include "Action.h"

class PolarizationVortexAction : public Action
{
public:
  PolarizationVortexAction(const std::string & name, InputParameters params);

  virtual void act();
private:
  const NonlinearVariableName _P_x;
  const NonlinearVariableName _P_y;
  const NonlinearVariableName _P_z;
  std::string _parameter;

};

template<>
InputParameters validParams<PolarizationVortexAction>();

#endif //POLARIZATIONVORTEXACTION_H
