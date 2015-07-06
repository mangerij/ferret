/**
 * @file   FerroelectricCouplingUWeak.h
 * @author J. Mangeri <mangerij@anl.gov>
 * @date   Jul. 1 2015
 * @brief  Implement the kernel for displacement variables corresponding to ferroelectic coupling energy,
 *         Assume the energy has the form -0.5*q_ijkl* ui_j * Pk_l where u is the displacement and P is the polarization.
 *         with zero jacobian
 */

#ifndef FERROELECTRICCOUPLINGUWEAK_H
#define FERROELECTRICCOUPLINGUWEAK_H

#include "Kernel.h"
#include "ElectrostrictiveTensorR4.h"

//Forward Declarations
class FerroelectricCouplingUWeak;

template<>
InputParameters validParams<FerroelectricCouplingUWeak>();

class FerroelectricCouplingUWeak: public Kernel
{
public:

  FerroelectricCouplingUWeak(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);



private:
  const MaterialProperty<ElectrostrictiveTensorR4> & _electrostrictive_tensor;
  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue& _polar_x;
  const VariableValue& _polar_y;
  const VariableValue& _polar_z;
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm

};
#endif //FERROELECTRICCOUPLINGUWEAK_H
