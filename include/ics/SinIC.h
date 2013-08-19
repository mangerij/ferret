/**
 * @file   ICTemplate.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Wed Jun 26 18:40:34 2013
 *
 * @brief
 *
 *
 */
#ifndef SINIC_H
#define SINIC_H

#include "Kernel.h"
#include "InitialCondition.h"
#include "InputParameters.h"

// Forward Declarations
class SinIC;

namespace libMesh { class Point; }

template<>
InputParameters validParams<SinIC>();


class SinIC:public InitialCondition
{
public:
  SinIC(const std::string & name, InputParameters parameters);
  virtual Real value(const Point & p);
private:
  Real _amplitude;
  Real _wave_length_x,_wave_length_y,_wave_length_z;
  Real _phrase_x,_phrase_y,_phrase_z;
  Real _vertical_shift;
};

#endif //SinIC
