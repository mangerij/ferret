/**
 * @file   SphereToCartFunc.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Aug 13 15:52:07 2013
 *
 * @brief
 *
 *
 */
#ifndef SPHERETOCARTFUNC_H
#define SPHERETOCARTFUNC_H

#include "Function.h"
#include "FunctionInterface.h"

class SphereToCartFunc;

template<>
InputParameters validParams<SphereToCartFunc>();

class SphereToCartFunc: public Function, FunctionInterface
{
public:
  SphereToCartFunc(const std::string & name, InputParameters parameters);

  virtual Real value(Real t, const Point & p);

protected:
  Function & _radial_func;
  Function & _polar_func;
  Function & _azimuthal_func;
  unsigned int _index;
};

#endif //SPHERETOCARTFUNC_H
