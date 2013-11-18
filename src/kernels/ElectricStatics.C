/**
 * @file   ElectricStatics.C
 * @author S. Gu <sgu@anl.gov>
 * @date   Tue Jun 11 10:08:26 2013
 *
 * @brief  Laplacian operator with permittivity
 *
 *
 */

#include "ElectricStatics.h"

template<>
InputParameters validParams<ElectricStatics>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("permittivity", "permittivity");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  params.addRequiredParam<Real>("polar_electric_scale","polar_electric scale");
  return params;
}



//Constructor
ElectricStatics::ElectricStatics(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _permittivity(getParam<Real>("permittivity")),
   _len_scale(getParam<Real>("len_scale")),
   _polar_electric_scale(getParam<Real>("polar_electric_scale"))
{
}


Real
ElectricStatics::computeQpResidual()
{
  return _permittivity*_grad_u[_qp]*_grad_test[_i][_qp]*_len_scale*_polar_electric_scale;
}

Real
ElectricStatics::computeQpJacobian()
{
   return _permittivity*_grad_phi[_j][_qp]*_grad_test[_i][_qp]*_len_scale*_polar_electric_scale;
}
