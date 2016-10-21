
#ifndef MATCHEDGRADVALUEBC_H
#define MATCHEDGRADVALUEBC_H

#include "IntegratedBC.h"

//Forward Declarations
class MatchedGradValueBC;

template<>
InputParameters validParams<MatchedGradValueBC>();

/**
 * Implements a simple coupled boundary condition where u=v on the boundary.
 */
class MatchedGradValueBC : public IntegratedBC
{
public:
  MatchedGradValueBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  const VariableGradient & _u_grad;
  const VariableValue & _v;

  /// The id of the coupled variable
  unsigned int _v_num;
  unsigned int _component;
};

#endif //MATCHEDVALUEBC_H
