#ifndef POLARIZATIONSURFACECHARGE_H
#define POLARIZATIONSURFACECHARGE_H

#include "IntegratedBC.h"
#include "FerretBase.h"

class PolarizationSurfaceCharge;

template<>
InputParameters validParams<PolarizationSurfaceCharge>();


class PolarizationSurfaceCharge : public FerretBase, public IntegratedBC
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  PolarizationSurfaceCharge(const InputParameters & parameters);


protected:
  virtual Real computeQpResidual();

  const VariableValue & _P_x;
  const VariableValue & _P_y;
  const VariableValue & _P_z;
  bool _J_polarization;
};


#endif //POLARIZATIONSURFACECHARGE_H
