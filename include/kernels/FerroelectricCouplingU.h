/**
 * @file   FerroelectricCouplingU.h
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 * @date   Jun. 15 2015
 * @brief  Implement the kernel for displacement variables corresponding to ferroelectic coupling energy,
 *         Assume the energy has the form -0.5*q_ijkl* ui_j * Pk_l where u is the displacement and P is the polarization.
 */

#ifndef FERROELECTRICCOUPLINGU_H
#define FERROELECTRICCOUPLINGU_H

#include "Kernel.h"
#include "ElectrostrictiveTensorR4.h"

//Forward Declarations
class FerroelectricCouplingU;

template<>
InputParameters validParams<FerroelectricCouplingU>();

class FerroelectricCouplingU: public Kernel
{
public:

  FerroelectricCouplingU(const std::string & name, InputParameters parameters);

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
#endif //FERROELECTRICCOUPLINGU_H
