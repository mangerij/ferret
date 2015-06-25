/**
 * @file   Electrostatics.h
 * @author S. Gu <sgu@anl.gov>
 * @modified J. Mangeri <mangerij@anl.gov>
 *
 */

#include "Electrostatics.h"

template<>
InputParameters validParams<Electrostatics>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("permittivity", "permittivity");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

//Constructor
Electrostatics::Electrostatics(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _permittivity(getParam<Real>("permittivity")),
   _len_scale(getParam<Real>("len_scale"))
{
}

Real
Electrostatics::computeQpResidual()
{
  return _permittivity * _grad_u[_qp] * _grad_test[_i][_qp] * _len_scale;
}

Real
Electrostatics::computeQpJacobian()
{
   return _permittivity * _grad_phi[_j][_qp] * _grad_test[_i][_qp] * _len_scale;
}
