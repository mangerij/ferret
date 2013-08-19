/**
 * @file   SinFunc.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Mon Aug 12 16:03:14 2013
 *
 * @brief compute sinusoidual wave function
 *
 *
 */
#ifndef SINFUNC_H
#define SINFUNC_H

#include "Function.h"

class SinFunc;

template<>
InputParameters validParams<SinFunc>();

class SinFunc: public Function
{
public:
  SinFunc(const std::string & name, InputParameters parameters);

  virtual Real value(Real t, const Point & p);

protected:
  Real _amplitude;
  Real _wave_length_x,_wave_length_y,_wave_length_z;
  Real _phrase_x,_phrase_y,_phrase_z;
  Real _vertical_shift;
};

#endif //SINFUNC_H
