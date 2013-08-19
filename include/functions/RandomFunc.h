/**
 * @file   RandomFunc.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Mon Aug 12 16:30:24 2013
 *
 * @brief
 *
 *
 */
#ifndef RANDOMFUNC_H
#define RANDOMFUNC_H

#include "Function.h"

class RandomFunc;

template<>
InputParameters validParams<RandomFunc>();

class RandomFunc: public Function
{
public:
  RandomFunc(const std::string & name, InputParameters parameters);

  virtual Real value(Real t, const Point & p);

protected:
  Real _min;
  Real _max;
  Real _range;
};

#endif //RANDOMFUNC_H
