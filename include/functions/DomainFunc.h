#ifndef DOMAINFUNC_H
#define DOMAINFUNC_H

/*
 Function initially distributes polarization
 to form c-domain structure.
 */

#include "Function.h"

class DomainFunc;

template <>
InputParameters validParams<DomainFunc>();

class DomainFunc : public Function
{
public:
  DomainFunc(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) override;

protected:
  Real _ax;
  Real _af;
  Real _min;
  Real _max;
};

#endif 
