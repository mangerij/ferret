/**
 * @file   ConstantLatticeMismatch.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 */

#ifndef CONSTANTLATTICEMISMATCH_H
#define CONSTANTLATTICEMISMATCH_H

#include "Kernel.h"

class ConstantLatticeMismatch;

template<>
InputParameters validParams<ConstantLatticeMismatch>();

class ConstantLatticeMismatch: public Kernel
{
public:

  ConstantLatticeMismatch(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const VariableValue & _Qxx;
  const VariableValue & _Qxy;
  const VariableValue & _Qxz;
  const VariableValue & _Qyx;
  const VariableValue & _Qyy;
  const VariableValue & _Qyz;
  const VariableValue & _Qzx;
  const VariableValue & _Qzy;
  const VariableValue & _Qzz;
  const unsigned int _disp_x_var;
  const unsigned int _disp_y_var;
  const unsigned int _disp_z_var;
  const VariableGradient & _disp_x_grad;
  const VariableGradient & _disp_y_grad;
  const VariableGradient & _disp_z_grad;
  const unsigned int _component;
  const unsigned int _deriv_component;
};
#endif
