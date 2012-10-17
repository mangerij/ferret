/****************************************************************/
/* Stress BC:                                                  */
/*     This BC is intended only for testing purpose.            */
/*                                                              */
/*     Anyway, it reads six values on designated boundary:      */
/*     stress_xx, stress_xy, stress_yy,                         */
/*     stress_yz, stress_zx, stress_zz                          */
/*     Then, it computes the traction by multiplying            */
/*     stress tensor with the normal to get the traction.       */
/*     The traction is used for the problem.                    */
/*                                                              */
/****************************************************************/

#ifndef STRESSBC_H
#define STRESSBC_H

#include "IntegratedBC.h"
//LibMesh includes
//#include "vector_value.h"

//Forward Declarations
class StressBC;

template<>
InputParameters validParams<StressBC>();

class StressBC : public IntegratedBC
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  StressBC(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();

  //private:
  const int _component; 
  const Real _stress_xx;
  const Real _stress_xy;
  const Real _stress_yy;
  const Real _stress_yz;
  const Real _stress_zx;
  const Real _stress_zz;
};

#endif
