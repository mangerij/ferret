/**
 * @file   WallEnergyDensity.h
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu Nov  7 12:03:22 2013
 *
 * @brief Calculate wall energy density
 *
 *
 */

#ifndef NORMALIZEDWALLENERGYDENSITY_H
#define NORMALIZEDWALLENERGYDENSITY_H

#include "AuxKernel.h"


//Forward Declarations
class NormalizedWallEnergyDensity;

template<>
InputParameters validParams<NormalizedWallEnergyDensity>();

/**
 * Coupled auxiliary value
 */
class NormalizedWallEnergyDensity : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  NormalizedWallEnergyDensity(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  const VariableGradient & _polar_x_grad;
  const VariableGradient & _polar_y_grad;
  const VariableGradient & _polar_z_grad;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const Real _G110,_G11, _G12, _G44, _G44P;
  const Real _len_scale;
};

#endif // NORMALIZEDWALLENERGYDENSITY_H
