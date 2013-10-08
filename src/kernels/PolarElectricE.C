/**
 * @file   PolarElectricE.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Thu May 30 11:59:56 2013
 *
 * @brief
 *
 *
 */

#include "PolarElectricE.h"

class PolarElectricE;

template<>
InputParameters validParams<PolarElectricE>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("permittivity", "permittivity");
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return params;
}



//Constructor
PolarElectricE::PolarElectricE(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _permittivity(getParam<Real>("permittivity")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z"))
{}



Real
PolarElectricE::computeQpResidual()
{
  return 0.5*(_polar_x[_qp]*_grad_test[_i][_qp](0)+_polar_y[_qp]*_grad_test[_i][_qp](1)+_polar_z[_qp]*_grad_test[_i][_qp](2));
}

Real
PolarElectricE::computeQpJacobian()
{
  return 0.0;
}

Real
PolarElectricE::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int component;
  if (jvar == _polar_x_var)
    component = 0;
  else if (jvar == _polar_y_var)
    component = 1;
  else if (jvar == _polar_z_var)
    component = 2;
  else
    return 0.0;
  return -0.5*_grad_test[_i][_qp](component)*_phi[_j][_qp];
}
