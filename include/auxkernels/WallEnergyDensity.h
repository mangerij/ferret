/**
 * @file   WallEnergyDensity.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu Nov  7 12:03:22 2013
 *
 * @brief Calculate wall energy density
 *
 *
 */

#ifndef WALLENERGYDENSITY_H
#define WALLENERGYDENSITY_H

#include "AuxKernel.h"


//Forward Declarations
class WallEnergyDensity;

template<>
InputParameters validParams<WallEnergyDensity>();

/**
 * Coupled auxiliary value
 */
class WallEnergyDensity : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  WallEnergyDensity(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  const VariableGradient & _polar_x_grad;
  const VariableGradient & _polar_y_grad;
  const VariableGradient & _polar_z_grad;
  const Real _G110,_G11, _G12, _G44, _G44P;
  const Real _len_scale;
};

#endif // WALLENERGYDENSITY_H
