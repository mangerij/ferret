/**
 * @file  ConstField.C
 * @author J. Mangeri <john.mangeri@uconn.edu>
 *
 * @brief test kernel for a constant field along z.
 */

#include "ConstField.h"

class ConstField;

template<>
InputParameters validParams<ConstField>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  params.addParam<Real>("len_scale", 1.0, "the len_scale of the unit");
  params.addParam<Real>("field", 0.0, "the constant field");
  return params;
}

ConstField::ConstField(const InputParameters & parameters)
  :Kernel(parameters),
   _polar_z_var(coupled("polar_z")),
   _polar_z(coupledValue("polar_z")),
   _len_scale(getParam<Real>("len_scale")),
   _field(getParam<Real>("field"))
{
}

Real
ConstField::computeQpResidual()
{
  return _field * _test[_i][_qp] * std::pow(_len_scale, 3.0);
}

Real
ConstField::computeQpJacobian()
{
  return 0.0;
}
