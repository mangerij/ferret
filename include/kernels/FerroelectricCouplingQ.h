/**
 * @file   FerroelectricCouplingQ.h
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov
 * @date   Jun 15. 2015

 * @brief  Implement the kernel for polar variables corresponding to ferroelectic coupling energy after
 *         the variational derivative of the polar dependent terms have been taken.
 *         This is only the quartic term. See notes.
 */
#ifndef FERROELECTRICCOUPLINGQ_H
#define FERROELECTRICCOUPLINGQ_H

#include "Kernel.h"
#include "ComputeElectrostrictiveTensor.h"
#include "Material.h"

//Forward Declarations
class FerroelectricCouplingQ;

template<>
InputParameters validParams<FerroelectricCouplingQ>();

class FerroelectricCouplingQ: public Kernel
{
public:

  FerroelectricCouplingQ(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const MaterialProperty<RankFourTensor> & _electrostrictive_tensor;
  const MaterialProperty<RankFourTensor> & _electrostrictive_tensorQ;
  const unsigned int _component;
  const unsigned int _polar_x_var;
  const unsigned int _polar_y_var;
  const unsigned int _polar_z_var;
  const VariableValue & _polar_x;
  const VariableValue & _polar_y;
  const VariableValue & _polar_z;
  const Real _artificial;
  const Real _len_scale;     //dimension unit, eg: 1e-9 for nm

};
#endif //FERROELECTRICCOUPLINGQ_H
