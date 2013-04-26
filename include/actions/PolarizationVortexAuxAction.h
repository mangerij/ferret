#ifndef POLARIZATIONVORTEXAUXACTION_H
#define POLARIZATIONVORTEXAUXACTION_H

#include "Action.h"
#include "FerretBase.h"

class PolarizationVortexAuxAction : public FerretBase, public Action
{
public:
  PolarizationVortexAuxAction(const std::string & name, InputParameters params);

  virtual void act();
private:
  const NonlinearVariableName _P_x;
  const NonlinearVariableName _P_y;
  const NonlinearVariableName _P_z;
  Real _a_x, _a_y, _c, _R, _L;
  std::string _p;
};

template<>
InputParameters validParams<PolarizationVortexAuxAction>();

#endif //POLARIZATIONVORTEXACTION_H
