/************************************************************************/
/* StressFree BC:                                                       */
/*     This BC is intended to implement sigma_i = 0                     */
/*     at a boundary for i = 0,1,2                                      */
/************************************************************************/

#ifndef STRESSFREEBC_H
#define STRESSFREEBC_H

#include "IntegratedBC.h"
#include "ComputeElectrostrictiveTensor.h"

//Forward Declarations
class StressFreeBC;

template<>
InputParameters validParams<StressFreeBC>();

class StressFreeBC : public IntegratedBC
{
public:

  StressFreeBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const MaterialProperty<RankFourTensor> & _electrostrictive_tensor;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  const unsigned int _component;
  const unsigned int _disp_x_var;
  const unsigned int _disp_y_var;
  const unsigned int _disp_z_var;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const VariableGradient & _disp_x_grad;
  const VariableGradient & _disp_y_grad;
  const VariableGradient & _disp_z_grad;
};

#endif /* STRESSFREEBC_H */
