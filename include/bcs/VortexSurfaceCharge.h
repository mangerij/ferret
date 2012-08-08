#ifndef VORTEX_SURFACE_CHARGE_H
#define VORTEX_SURFACE_CHARGE_H

#include "IntegratedBC.h"


class VortexSurfaceCharge;

template<>
InputParameters validParams<VortexSurfaceCharge>();

/**
 * Implements a Neumann boundary condition of electrostatics corresponding
 * to a planar polarization vortex located in a cylindrical disk (dot)
 * at distance a along the x-axis from the dot center. This is basically 
 * NeumannBC with a position-dependent value, assuming the boundary is
 * the lateral (vertical) surface of a cylinder (that, the normal is 
 * horizontal, as having small 'verticality' -- its z-component). 
 */
class VortexSurfaceCharge : public IntegratedBC
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  VortexSurfaceCharge(const std::string & name, InputParameters parameters);


protected:
  virtual Real computeQpResidual();

  Real _a;
  Real _verticality;
};


#endif //VORTEX_SURFACE_CHARGE_H
