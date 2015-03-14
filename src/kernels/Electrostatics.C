/**
 * @file   ElectricStatics.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jun 11 10:08:26 2013
 *
 * @brief  Laplacian operator with permittivity
 *
 *
 */

#include "Electrostatics.h"

template<>
InputParameters validParams<Electrostatics>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("permittivity", "permittivity");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  params.addParam<Real>("polar_electric_scale",1.0,"polar_electric scale");
  return params;
}



//Constructor
Electrostatics::Electrostatics(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _permittivity(getParam<Real>("permittivity")),
   _len_scale(getParam<Real>("len_scale")),
   _polar_electric_scale(getParam<Real>("polar_electric_scale"))
{
}


Real
Electrostatics::computeQpResidual()
{
  return _permittivity*_grad_u[_qp]*_grad_test[_i][_qp]*_len_scale*_polar_electric_scale;
}

Real
Electrostatics::computeQpJacobian()
{
   return _permittivity*_grad_phi[_j][_qp]*_grad_test[_i][_qp]*_len_scale*_polar_electric_scale;
}
