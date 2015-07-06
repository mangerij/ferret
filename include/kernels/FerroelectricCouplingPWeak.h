/**
 * @file   FerroelectricCouplingPWeak.h
 * @author J. Mangeri <mangerij@anl.gov
 * @date   Jul 1. 2015
 * @brief   Implement the kernel for polar variables corresponding to ferroelectic coupling energy.
 *         Assume the energy has the form -0.5*q_ijkl* ui_j * Pk_l where u is the displacement and P is the polarization.
 *         with a zero jacobian
 */
#ifndef FERROELECTRICCOUPLINGPWEAK_H
#define FERROELECTRICCOUPLINGPWEAK_H

#include "Kernel.h"
#include "ElectrostrictiveTensorR4.h"

//Forward Declarations
class FerroelectricCouplingPWeak;

template<>
InputParameters validParams<FerroelectricCouplingPWeak>();

class FerroelectricCouplingPWeak: public Kernel
{
public:

  FerroelectricCouplingPWeak(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);


private:
  const MaterialProperty<ElectrostrictiveTensorR4> & _electrostrictive_tensor;
  const unsigned int _component;
  const unsigned int _disp_x_var;
  const unsigned int _disp_y_var;
  const unsigned int _disp_z_var;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableGradient& _disp_x_grad;
  const VariableGradient& _disp_y_grad;
  const VariableGradient& _disp_z_grad;
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm

};
#endif //FERROELECTRICCOUPLINGPWEAK_H
