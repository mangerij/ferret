/****************************************************************/
/* Hydrostatic Dirichlet BC:                                    */
/*     This BC is for a sphere....                              */
/*                                                              */
/****************************************************************/

#ifndef HYDROSTATICDIRICHLETBC_H
#define HYDROSTATICDIRICHLETBC_H

#include "NodalBC.h"
//LibMesh includes
//#include "libmesh/vector_value.h"

//Forward Declarations
class HydrostaticDirichletBC;

template<>
InputParameters validParams<HydrostaticDirichletBC>();

class HydrostaticDirichletBC : public NodalBC
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  HydrostaticDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  //private:
  const Real _center_x;
  const Real _center_y;
  const Real _center_z;
  const Real _displacement;
  const int  _component;
};

#endif
