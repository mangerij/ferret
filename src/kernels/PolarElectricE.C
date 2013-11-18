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
  params.addRequiredParam<Real>("polar_electric_scale","polar_electric scale");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");

  return params;
}



//Constructor
PolarElectricE::PolarElectricE(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _permittivity(getParam<Real>("permittivity")),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValueOld("polar_x")),
   _polar_y(coupledValueOld("polar_y")),
   _polar_z(coupledValueOld("polar_z")),
   _len_scale(getParam<Real>("len_scale")),
   _polar_electric_scale(getParam<Real>("polar_electric_scale"))
{}



Real
PolarElectricE::computeQpResidual()
{
  return -((_polar_x[_qp]*_grad_test[_i][_qp](0)+_polar_y[_qp]*_grad_test[_i][_qp](1)+_polar_z[_qp]*_grad_test[_i][_qp](2))*pow(_len_scale,2.0))*_polar_electric_scale;
}

Real
PolarElectricE::computeQpJacobian()
{
  return 0.0;
}

Real
PolarElectricE::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
