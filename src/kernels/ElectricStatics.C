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
  return params;
}



//Constructor
ElectricStatics::ElectricStatics(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _permittivity(getParam<Real>("permittivity"))
{
}


Real
ElectricStatics::computeQpResidual()
{
  return _permittivity*_grad_u[_qp]*_grad_test[_i][_qp];
}

Real
ElectricStatics::computeQpJacobian()
{
   return _permittivity*_grad_phi[_j][_qp]*_grad_test[_i][_qp]; 
}
