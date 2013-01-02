/****************************************************************/
/* Hydrostatics BC:                                             */
/*     This BC applies a hydrostatic pressur to  sphere.        */
/*     The pressure is specified.                               */ 
/*     Then, the traction is computed by taking the diagonals of*/
/*     stress tensor to be a vector parallel to the surface normal. */
/*     The traction is obtained by multiplying the              */
/*     stress tensor with the normal to get the traction.       */
/*     The traction is used for the problem.                    */
/*                                                              */
/****************************************************************/

#ifndef HYDROSTATICBC_H
#define HYDROSTATICBC_H

#include "IntegratedBC.h"

//Forward Declarations
class HydrostaticBC;

template<>
InputParameters validParams<HydrostaticBC>();

class HydrostaticBC : public IntegratedBC
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  HydrostaticBC(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();

private:
  Real _pressure;
  int _component;
};

#endif
