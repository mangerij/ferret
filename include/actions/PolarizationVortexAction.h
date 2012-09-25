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
  Real _a_x, _a_y, _c, _R, _L;
  std::string _p;
  bool _debug;
};

template<>
InputParameters validParams<PolarizationVortexAction>();

#endif //POLARIZATIONVORTEXACTION_H
