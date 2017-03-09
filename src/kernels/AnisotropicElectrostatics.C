/**
 * @file   AnisotropicElectrostatics.h
 * @author J. Mangeri <john.mangeri@uconn.edu>
 */

#include "AnisotropicElectrostatics.h"

template<>
InputParameters validParams<AnisotropicElectrostatics>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("inplane_permittivity", "in plane permittivity");
  params.addRequiredParam<Real>("outofplane_permittivity", "out of plane permittivity");
  params.addParam<Real>("len_scale", 1.0, "the length scale of the unit");
  return params;
}

AnisotropicElectrostatics::AnisotropicElectrostatics(const InputParameters & parameters)
  :Kernel(parameters),
   _inplane_permittivity(getParam<Real>("inplane_permittivity")),
   _outofplane_permittivity(getParam<Real>("outofplane_permittivity")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
AnisotropicElectrostatics::computeQpResidual()
{
  return _inplane_permittivity * (_grad_u[_qp](0) * _grad_test[_i][_qp](0) + (_outofplane_permittivity/_inplane_permittivity) * _grad_u[_qp](1) * _grad_test[_i][_qp](1) + _grad_u[_qp](2) * _grad_test[_i][_qp](2) ) * _len_scale;
}

Real
AnisotropicElectrostatics::computeQpJacobian()
{
  return _inplane_permittivity * (_grad_phi[_j][_qp](0) * _grad_test[_i][_qp](0) + (_outofplane_permittivity/_inplane_permittivity) * _grad_phi[_j][_qp](1) * _grad_test[_i][_qp](1) + _grad_phi[_j][_qp](2) * _grad_test[_i][_qp](2) ) * _len_scale;
}
