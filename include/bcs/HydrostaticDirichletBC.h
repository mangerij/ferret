/****************************************************************/
/* Hydrostatic Dirichlet BC:                                    */
/*     This BC is for a sphere....                              */
/*                                                              */
/****************************************************************/

#ifndef HYDROSTATICDIRICHLETBC_H
#define HYDROSTATICDIRICHLETBC_H

#include "NodalBC.h"
//LibMesh includes
//#include "vector_value.h"

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
  HydrostaticDirichletBC(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();

  //private:
  const int  _component;
  const Real _displacement; 
  const Real _center_x;
  const Real _center_y;
  const Real _center_z;
};

#endif
