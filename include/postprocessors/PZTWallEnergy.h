/**
* @file   PZTWallEnergy.C
* @author J. Mangeri <john.mangeri@uconn.edu>
* @date   Dec 15 2015
 *
 * @brief
 *
 *
 */

#ifndef PZTWALLENERGY_H
#define PZTWALLENERGY_H

#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class PZTWallEnergy;

template<>
InputParameters validParams<PZTWallEnergy>();

class PZTWallEnergy : public ElementIntegralPostprocessor
{
public:
  PZTWallEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableGradient& _polar_x_grad;
  const VariableGradient& _polar_y_grad;
  const VariableGradient& _polar_z_grad;
  const Real _G110;
  const Real _len_scale;
};

#endif
