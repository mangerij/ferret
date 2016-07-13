/**
 * @file   PolarElectricEStrong.C
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 *
 */

#include "PolarElectricEStrong.h"

class PolarElectricEStrong;

template<>
InputParameters validParams<PolarElectricEStrong>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  params.addCoupledVar("polar_y", 0.0, "The y component of the polarization");
  params.addCoupledVar("polar_z", 0.0, "The z component of the polarization");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

PolarElectricEStrong::PolarElectricEStrong(const InputParameters & parameters)
  :Kernel(parameters),
   _polar_x_var(coupled("polar_x")),
   _polar_y_var(coupled("polar_y")),
   _polar_z_var(coupled("polar_z")),
   _polar_x(coupledValue("polar_x")),
   _polar_y(coupledValue("polar_y")),
   _polar_z(coupledValue("polar_z")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
PolarElectricEStrong::computeQpResidual()
{
  Real RpolarE = 0.0;
  RpolarE += - (_polar_x[_qp] * _grad_test[_i][_qp](0) + _polar_y[_qp] * _grad_test[_i][_qp](1) + _polar_z[_qp] * _grad_test[_i][_qp](2)) * std::pow(_len_scale, 2.0);
  ///  Moose::out << "\n R_polarE-"; std::cout << " = " << RpolarE;
  return RpolarE;
}
Real
PolarElectricEStrong::computeQpJacobian()
{
  return 0.0;
}

Real
PolarElectricEStrong::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _polar_x_var)
    return - _phi[_j][_qp] * _grad_test[_i][_qp](0) * std::pow(_len_scale, 2.0);
  else if (jvar == _polar_y_var)
    return - _phi[_j][_qp] * _grad_test[_i][_qp](1) * std::pow(_len_scale, 2.0);
  else if (jvar == _polar_z_var)
    return - _phi[_j][_qp] * _grad_test[_i][_qp](2) * std::pow(_len_scale, 2.0);
  else
    return 0.0;
}
