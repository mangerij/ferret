/**
 * @file   SinFuncAux.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Mon Aug 12 12:12:20 2013
 *
 * @brief explicitly compute sinusoidual wave function
 *
 *
 */

#ifndef SINFUNCAUX_H
#define SINFUNCAUX_H

#include "AuxKernel.h"


//Forward Declarations
class SinFuncAux;

template<>
InputParameters validParams<SinFuncAux>();

/**
 * Coupled auxiliary value
 */
class SinFuncAux : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  SinFuncAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();
  Real _amplitude;
  Real _wave_length_x,_wave_length_y,_wave_length_z;
  Real _phrase_x,_phrase_y,_phrase_z;
  Real _vertical_shift;
};

#endif // SINFUNCAUX_H
