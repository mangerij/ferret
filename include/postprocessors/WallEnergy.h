/**
 * @file   WallEnergy.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jul 30 13:09:54 2013
 *
 * @brief
 *
 *
 */

#ifndef WALLENERGY_H
#define WALLENERGY_H

//TODO: include the base header
#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class WallEnergy;

template<>
InputParameters validParams<WallEnergy>();

//TODO: change the base class!
class WallEnergy : public ElementIntegralPostprocessor
{
public:
  WallEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  const VariableGradient& _polar_x_grad;
  const VariableGradient& _polar_y_grad;
  const VariableGradient& _polar_z_grad;
  const Real _G110,_G11, _G12, _G44, _G44P;
  const Real _len_scale;
};

#endif
