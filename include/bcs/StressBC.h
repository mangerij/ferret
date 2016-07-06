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
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

//LibMesh includes
//#include "libmesh/vector_value.h"

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
  StressBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  std::vector<Real> _stress_vector;

  RankTwoTensor _boundary_stress;
  std::vector<const VariableValue *> _boundary_stress_vars;

  const MaterialProperty<RankFourTensor> & _Jacobian_mult;


  const int _component;

  bool _convert_to_gpa;

  Real _multiplier;
};

#endif
