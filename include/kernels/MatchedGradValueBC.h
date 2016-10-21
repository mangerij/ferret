#ifndef MATCHEDVALUEBC_H
#define MATCHEDVALUEBC_H

#include "NodalBC.h"

//Forward Declarations
class MatchedGradValueBC;

template<>
InputParameters validParams<MatchedGradValueBC>();

/**
 * Implements a simple coupled boundary condition where u=v on the boundary.
 */
class MatchedGradValueBC : public NodalBC
{
public:
  MatchedGradValueBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _v;

  /// The id of the coupled variable
  unsigned int _v_num;
};

#endif //MATCHEDGRADVALUEBC_H
