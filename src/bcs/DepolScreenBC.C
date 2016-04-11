/* Implement a simple screening of the depolarizing field.
   More complicated methods will come later. */

#include "DepolScreenBC.h"

template<>
InputParameters validParams<DepolScreenBC>()
{
  InputParameters p = validParams<NodalNormalBC>();
  p.addRequiredParam<Real>("value", "Value of the BC");
  p.addRequiredCoupledVar("polar_x", "The x component of the polarization");
  p.addRequiredCoupledVar("polar_y", "The y component of the polarization");
  p.addRequiredCoupledVar("polar_z", "The z component of the polarization");
  return p;
}

DepolScreenBC::DepolScreenBC(const InputParameters & parameters) :
  NodalNormalBC(parameters),
  _polar_x(coupledValue("polar_x")),
  _polar_y(coupledValue("polar_y")),
  _polar_z(coupledValue("polar_z")),
  _value(getParam<Real>("value"))
{}

Real
DepolScreenBC::computeQpResidual()
{
  return _polar_x[_qp] * _normal(0) + _polar_y[_qp] * _normal(1) + _polar_z[_qp] * _normal(2) - _value * (_polar_x[_qp] * _normal(0) + _polar_y[_qp] * _normal(1) + _polar_z[_qp] * _normal(2));
}
