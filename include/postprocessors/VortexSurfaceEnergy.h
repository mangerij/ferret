#ifndef VORTEX_SURFACE_ENERGY_H
#define VORTEX_SURFACE_ENERGY_H

#include "SideIntegral.h"

//Forward Declarations
class VortexSurfaceEnergy;

template<>
InputParameters validParams<VortexSurfaceEnergy>();

/**
 * This postprocessor computes a surface integral of the product of the potential against 
 * a vortex charge.
 */
class VortexSurfaceEnergy : public SideIntegral
{
public:
  VortexSurfaceEnergy(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpIntegral();
private:
  Real _a;
  Real _verticality;
};

#endif
